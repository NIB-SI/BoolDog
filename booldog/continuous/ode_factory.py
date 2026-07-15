'''Build ODE systems (continuous relaxations) from Boolean networks.

:func:`ode_factory` is the main entry point: given a
:py:class:`booldog.BoolDogModel` and a transform name, it returns an
:py:class:`ODE` subclass instance (:py:class:`BooleCubeODE` or
:py:class:`SquadODE`) whose `dxdt` callable is suitable for
`scipy.integrate.solve_ivp`.
'''

import inspect
import logging
from itertools import product
import numpy as np


from booldog.utils import parameter_to_array, ensure_ndarray

logger = logging.getLogger(__name__)

##############################
#      GENERAL FUNCTIONS     #
##############################


def ode_factory(network, transform, **kwargs):
    '''Create an :py:class:`ODE` instance from a Boolean network.

    Parameters
    ----------
    network : booldog.BoolDogModel
        Input Boolean network to transform.

    transform : str
        One of the accepted transforms (case-insensitive). See
        :py:data:`transforms` or :ref:`Notes <tagnotesode>` for options.

    Other Parameters
    ----------------
    **kwargs
        Additional arguments and keyword arguments passed to the
        selected ODE class's constructor.

    Returns
    -------
    ode : ODE
        A :py:class:`BooleCubeODE` or :py:class:`SquadODE` instance
        (an :py:class:`ODE` subclass), depending on `transform`.


    .. _tagnotesode:

    Notes
    -----
    For specific transforms, see the relevant class for keyword
    arguments (`**kwargs`).

    The class per transform is defined in
    :py:data:`ode_classes`.

    If the parameter is an int or float, the value is assigned for
    all variables. Otherwise the parameter argument should be a dict
    with keys as node names and values for their initial state. In
    this case, if the initial state is not defined for all nodes, a
    `default` key with the default value should also be present in the
    dict.

    Here follows a summary of the transform-specific keyword arguments

        'squad'

            - :py:class:`SquadODE`
            - gamma : decay rate
            - h : sigmoidal gain

        'boolecube'

            - :py:class:`BooleCubeODE`
            - tau : life-time of species

        'hillcube'

            - :py:class:`BooleCubeODE`
            - tau : life-time of species
            - n : Hill coefficient
            - k : Hill dissociation constant

        'normalisedhillcube'

            - :py:class:`BooleCubeODE`
            - tau : life-time of species
            - n : Hill coefficient
            - k : Hill dissociation constant
    '''
    transform = transform.lower()
    if transform == 'placeholder':
        class_ = BooleCubeODE
    elif not transform in ode_classes:
        raise ValueError(f"transform' argument must be one of"\
                         f"{list(ode_classes.keys())}")
    else:
        class_ = ode_classes[transform]

    # for sphinx documentation
    # ODE_factory.ex_class = ODE
    # ODE_factory.ex_class.__bases__ = tuple(set(ode_classes.values()))

    return class_(network, transform, **kwargs)


##############################
#        CHILD CLASS         #
##############################


class ODE():
    '''Base class for the ODE systems built by :func:`ode_factory`.

    Not intended to be instantiated directly. Its subclasses
    (:py:class:`BooleCubeODE`, :py:class:`SquadODE` — selected via the
    `transform` argument) each build their own `dxdt` right-hand-side
    function; this class only provides the shared `event_function` (used
    by `scipy.integrate.solve_ivp` to stop integration exactly at a
    perturbation time) and `update` (rebuild `dxdt`, e.g. after freezing
    a subset of nodes) machinery.
    '''

    def __init__(self, network, transform):
        '''Initialise the base ODE attributes.

        Parameters
        ----------
        network : booldog.BoolDogModel
            Input Boolean network.
        transform : str
            Name of the transform being constructed (e.g. ``'squad'``,
            ``'boolecube'``); stored as-is on `self.transform`.

        Notes
        -----
        If `transform` is ``'placeholder'`` this returns immediately
        without setting any attributes.
        '''
        if transform == 'placeholder':
            return

        # if not isinstance(network, RegulatoryNetwork):
        #     raise TypeError(f"'network' argument must be a RegulatoryNetwork object."\
        #                     f"not {type(network)}. ")

        self.n = len(network)
        '''int : number of nodes in the network.'''
        self.boolean_network = network
        '''booldog.BoolDogModel : the underlying Boolean network.'''
        self.transform = transform
        '''str : name of the transform used to build this ODE system.'''

        logger.info("Creating ODE system for %s.", transform)

    def event_function(self, t, x, event_t, *args):
        '''Event function for `events` of `scipy.integrate.solve_ivp`.

        Returns the signed time-to-event (`t - event_t`); its root is
        `event_t`, so `solve_ivp` uses it to stop integration exactly at
        that time (e.g. so a node perturbation can then be applied).

        Parameters
        ----------
        t : float
            Current time-point of the simulation.
        x : ndarray
            Current state (unused; present to match the event-function
            signature expected by `solve_ivp`).
        event_t : float
            Time-point at which the event should trigger.
        *args
            Ignored; absorbs any extra positional arguments `solve_ivp`
            passes via its `args=` parameter.
        '''
        return t - event_t

    event_function.terminal = True
    '''bool : Marks `event_function` as a terminal event, i.e.
    `scipy.integrate.solve_ivp` stops integration as soon as it
    triggers.'''

    def update(self, off_nodes=None):
        ''' Recompute and reset `self.dxdt`.

        Shortcut to the ODE subclass's `_get_system` method; used e.g.
        to change which nodes are held constant partway through a
        simulation.

        Parameters
        ----------
        off_nodes : list of int, optional
            List of node **indices** to set derivative to 0,
            i.e. these nodes will remain constant in simulation.
        '''
        self.dxdt = self._get_system(off_nodes=off_nodes)


##############################
#       PARENT CLASSES       #
##############################


# https://github.com/krumsieklab/Odefy/blob/11d048d550a8f64250ba01f76f5a83048c8be6cf/Odefy-1.20/code/models/CreateCubeCalls.m
class BooleCubeODE(ODE):
    '''ODE system built via multivariate polynomial (multilinear)
    interpolation of a Boolean network's prime implicants.

    Backs the ``'boolecube'``, ``'hillcube'``, and
    ``'normalisedhillcube'`` transforms; `transform` selects which
    `transform_function` (`identity`, `hill`, or `normalised_hill`
    respectively) is applied to the state before interpolation.
    '''

    def __init__(self, network, transform, tau=1, n=3, k=0.5, **kwargs):
        ''' Initialise a BooleCubeODE system.

        Parameters
        ----------
        network : booldog.BoolDogModel
            Input Boolean network.
        transform : str
            One of ``'boolecube'``, ``'hillcube'``,
            ``'normalisedhillcube'``; selects `transform_function`.

        tau : int, float, or dict, optional
            Life-time of species (per-node relaxation time). A node
            with ``tau == 0`` has its derivative forced to 0, i.e. it is
            held constant.

        n : int, float, or dict, optional
            Hill coefficient, used by the `hill`/`normalised_hill`
            transform functions (has no effect for ``'boolecube'``,
            which uses `identity`).

        k : int, float, or dict, optional
            Hill dissociation constant, used by the `hill`/
            `normalised_hill` transform functions (has no effect for
            ``'boolecube'``).

        Other Parameters
        ----------------
        **kwargs
            Ignored; accepted so that :func:`ode_factory` can pass
            shared keyword arguments without every transform needing to
            accept them.

        References
        ----------
        [1] Wittmann, D. M., Krumsiek, J., Saez-Rodriguez, J., Lauffenburger,
        D., A., Klamt, S., & Theis, F. J. (2009). Transforming Boolean models
        to continuous models: Methodology and application to T-cell receptor
        signaling. BMC Systems Biology, 3(1), 98.
        https://doi.org/10.1186/1752-0509-3-98

        '''
        super().__init__(network, transform)


        self.param_n = parameter_to_array(n, self.boolean_network.index)
        '''arraylike : Hill coefficient'''
        self.param_k = parameter_to_array(k, self.boolean_network.index)
        '''arraylike : Hill dissociation constant'''

        self.param_tau = parameter_to_array(tau, self.boolean_network.index)
        '''arraylike : life-time of species'''

        self.param_dict = {
            "n": self.param_n,
            "k": self.param_k,
            "tau": self.param_tau
        }
        '''dict : convenience mapping from parameter name ("n", "k",
        "tau") to its per-node array (`param_n`, `param_k`,
        `param_tau`).'''

        if transform == 'boolecube':
            self.transform_function = self.identity

        elif transform == 'hillcube':
            self.transform_function = self.hill

        elif transform == 'normalisedhillcube':
            self.transform_function = self.normalised_hill

        else:
            raise TypeError(f"Unknown transform {transform}. ")

        # returns an array function
        self.B1 = self.homologue_b1()
        '''function : multilinear interpolation of the network's Boolean
        rules, ``B1(x) -> ndarray``; see `homologue_b1`.'''

        # returns an array function
        self.dxdt = self._get_system()
        '''function : the ODE system's right-hand side, dx/dt = f(t, x)'''


    def hill(self, x_array):
        '''Hill-function transform of a state vector.

        Parameters
        ----------
        x_array : ndarray
            State to transform, shape ``(n,)`` (per-node values, typically
            in [0, 1]).

        Returns
        -------
        ndarray
            ``x_array**n / (x_array**n + k**n)``, computed element-wise
            using the per-node `param_n`/`param_k` arrays. Same shape as
            `x_array`.
        '''
        return x_array**self.param_n / \
               (x_array**self.param_n + self.param_k**self.param_n)

    def normalised_hill(self, x_array):
        '''Hill-function transform, normalised so that
        ``normalised_hill(1) == 1``.

        Parameters
        ----------
        x_array : ndarray
            State to transform, shape ``(n,)``.

        Returns
        -------
        ndarray
            `hill(x_array) / hill(1)`, element-wise. Same shape as `x_array`.
        '''
        return self.hill(x_array) / self.hill(1)

    def identity(self, x_array):
        '''No-op transform, used for the ``'boolecube'`` transform.

        Parameters
        ----------
        x_array : ndarray
            State to transform, shape ``(n,)``.

        Returns
        -------
        ndarray
            `x_array`, unchanged.
        '''
        return x_array

    def _get_system(self, off_nodes=None):
        '''Build the `dxdt(t, x_array, *args)` right-hand-side function.

        Nodes with `param_tau == 0` are automatically added to
        `off_nodes` (their derivative is always 0), in addition to any
        indices passed in explicitly.

        Parameters
        ----------
        off_nodes : iterable of int, optional
            Node **indices** whose derivative should be forced to 0
            (held constant), on top of any zero-`tau` nodes.

        Returns
        -------
        dxdt : callable
            Function ``dxdt(t, x_array, *args) -> ndarray`` giving
            ``dx/dt = (B1(transform_function(x)) - x) / tau``
            element-wise, after clipping `x_array` to [0, 1]; suitable
            as the `fun` argument of `scipy.integrate.solve_ivp`.
        '''

        if off_nodes is None:
            off_nodes = set()
        else:
            off_nodes = set(off_nodes)

        off_nodes.update(np.where(self.param_tau == 0)[0])

        def dxdt(t, x_array, *args):
            x_array[x_array < 0] = 0
            x_array[x_array > 1] = 1

            b = self.B1(self.transform_function(x_array))
            d = 1 / self.param_tau * (b-x_array)
            for i in off_nodes:
                d[i] = 0
            return d

        return dxdt

    def homologue_b1(self):
        '''Build the multilinear-interpolation function `B1(x)`.

        For each node, enumerates the states consistent with its
        positive prime implicants (expanding any parents left free by a
        prime over all their combinations) and sums, over those states,
        a product of ``x[j]`` (parent j on in that state) or
        ``(1-x[j])`` (parent j off) terms - the standard multilinear
        extension of a Boolean function to the unit hypercube [1]. The
        resulting per-node expression strings are compiled with `eval`
        into a single vectorised function operating on the whole state
        array.

        If a node's expression is too large to compile (`eval` raising
        `RecursionError`), that node's contribution is left as the
        constant string ``'0'`` and a message is logged.

        Returns
        -------
        B1 : callable
            Function ``B1(x) -> ndarray`` mapping a state vector ``x``
            (length ``self.n``) to the per-node multilinear-interpolation
            values.

        References
        ----------
        [1] Wittmann et al. (2009); see the references on
        `BooleCubeODE`.
        '''
        # spaces = set()
        # sums = []
        # all_B1s = ['0']*self.boolean_network.n
        # for node in self.boolean_network.nodes: # iterate over all nodes
        #     for prime_dict in self.boolean_network.primes[node][1]:
        #         for x_bool in self.boolean_network.generate_states(
        #                                     fixed=prime_dict):
        #             if not (tuple(x_bool) in spaces):
        #                 spaces.add(tuple(x_bool))
        #                 product = []
        #                 for i, b in enumerate(x_bool):
        #                     if b ==0:
        #                         product.append(f'(1-x[{i}])')
        #                     else:
        #                         product.append(f'x[{i}]')
        #             sums.append('*'.join(product))
        #     B1 = " + ".join(sums)
        #     if B1 != '':
        #         all_B1s[self.boolean_network.index[node]]  = B1

        all_B1s = ['0'] * self.boolean_network.n
        for node in self.boolean_network.nodes:  # iterate over all nodes
            parents = set(self.boolean_network.get_parents(node))

            states = []
            spaces = set()
            for prime_dict in self.boolean_network.primes[node][1]:

                fixed_parents = prime_dict.keys()
                free_parents = parents - set(fixed_parents)

                if len(free_parents) > 0:
                    for x in product([0, 1], repeat=len(free_parents)):
                        this_state = prime_dict.copy()
                        for parent, parent_state in zip(free_parents, x):
                            this_state[parent] = parent_state
                        str_rep = "".join([
                            str(this_state[node])
                            if node in this_state else "-"
                            for node in self.boolean_network.nodes
                        ])
                        if not str_rep in spaces:
                            states.append(this_state)
                            spaces.add(str_rep)
                else:
                    states.append(prime_dict)

            sums = []
            str_sums = []

            terms = []
            str_terms = []

            for state in states:
                subterms = []
                str_subterms = []
                for this_node, this_node_state in state.items():
                    if this_node_state == 0:
                        subterms.append(
                            f'(1-x[{self.boolean_network.index[this_node]}])')
                        str_subterms.append(f'(1-x[{this_node}])')
                    else:
                        subterms.append(
                            f'x[{self.boolean_network.index[this_node]}]')
                        str_subterms.append(f'   x[{this_node}] ')
                term = '*'.join(subterms)
                str_term = ' * '.join(str_subterms)

                sums.append(term)
                str_sums.append(str_term)

            B1 = " + ".join(sums)
            str_B1 = "\n + ".join(str_sums)

            if B1 != '':
                try:
                    eval(f'lambda x:{B1}')
                except RecursionError:
                    B1 = ''
                    logger.info(
                        "The rule for node %s is too long (depends on %i states. ",
                        node, len(states))

                all_B1s[self.boolean_network.index[node]] = B1

            logger.debug('%s %i', node, len(states))

        return eval('lambda x:' + 'np.array([' + ','.join(all_B1s) + '])')
        #return lambda x: x

    def write_c_code(self):
        '''Not implemented.

        Raises
        ------
        NotImplementedError
        '''
        raise NotImplementedError("TODO!")


class SquadODE(ODE):
    '''ODE system built via the SQUAD sigmoidal transform of a Boolean
    network's activator/inhibitor structure.

    Backs the ``'squad'`` transform.

    Attributes
    ----------
    activations : arraylike
        n x n activator matrix (`activations[i, j] == 1` iff node j
        activates node i); see `booldog.boolean.boolean.
        BooleanNetworkMixin.primes_to_matrices`.

    inhibitions : arraylike
        n x n inhibitor matrix (`inhibitions[i, j] == 1` iff node j
        inhibits node i); see `primes_to_matrices`.

    '''

    def __init__(self, network, transform, gamma=1, h=10, **kwargs):
        '''Build a SQUAD ODE system from a Boolean network's activator/
        inhibitor structure.

        Parameters
        ----------
        network : booldog.BoolDogModel
            Input Boolean network.

        transform : str
            Expected to be ``'squad'``.

        gamma : int, float, or dict, optional
            Per-node decay rate.

        h : int, float, or dict, optional
            Per-node sigmoidal gain.

        Other Parameters
        ----------------
        **kwargs
            Ignored; accepted so that :func:`ode_factory` can pass
            shared keyword arguments without every transform needing to
            accept them.

        References
        ----------
        [1] Di Cara, A., Garg, A., De Micheli, G., Xenarios, I., & Mendoza, L.
        (2007). Dynamic simulation of regulatory networks using SQUAD.
        BMC Bioinformatics, 8(1), 1–10. https://doi.org/10.1186/1471-2105-8-462
        '''
        super().__init__(network, transform)

        self.param_gamma = parameter_to_array(gamma,
                                              self.boolean_network.index)
        '''arraylike : decay rate'''
        self.param_h = parameter_to_array(h, self.boolean_network.index)
        '''arraylike : sigmoidal gain'''

        # print(self.param_gamma)
        # print(self.param_h)
        # matrices
        self.activations, self.inhibitions = self.boolean_network.primes_to_matrices()

        # needed for computations
        col_ones = np.ones((self.n))
        self._A1 = self.activations.dot(col_ones)
        # Di Cara et al (2007) eq. 2 has (1 - sum(alpha_n)) / sum(alpha_n)
        # for the activator prefactor (a "-", not a "+"). Probably an
        # error: taken literally it violates the paper's own 0 <= omega <= 1
        # bound (e.g. two fully-active default-weight activators give omega =
        # -1/3), whereas "+1" here (as implemented below) correctly gives
        # omega = 1 when fully activated with no inhibitors.
        self._a1 = (1 + self._A1) / self._A1
        self._B1 = self.inhibitions.dot(col_ones)
        self._b1 = (1 + self._B1) / self._B1

        self.dxdt = self._get_system()
        '''function : the ODE system's right-hand side, dx/dt = f(t, x)'''

    def _omega(self, x):
        '''
        Equation (2) of Di Cara et al (2007)
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-462
        Based on Andre Blejec R code.
        '''

        x = ensure_ndarray(x)

        a_x = self.activations.dot(x)
        a = ensure_ndarray(self._a1 * a_x / (1+a_x))
        a[~np.isfinite(a)] = 1

        b_x = self.inhibitions.dot(x)
        b = ensure_ndarray(self._b1 * b_x / (1+b_x))
        b[~np.isfinite(b)] = 0

        o = ensure_ndarray(a * (1-b))
        o[np.where(self._A1 + self._B1 == 0)] = 0

        return o

    def _dxdt_transform(self, x, w):
        ''' Equation (2) of Di Cara et al (2007).

        The numerator's second exponent here uses (w - 0.5); the 2007 paper's
        equation has just w (no shift). This matches the boundary-normalised
        variant used in later work building on SQUAD (e.g. Martinez-Sosa &
        Mendoza, 2013), which gives f(0)=0 and f(1)=1 exactly whilethe 2007
        form does not (f(0) != 0, f(1) != 1 for finite h).
        Confirmed against a later SQUAD reimplementation in
        https://github.com/caramirezal/SQUADBookChapter
        '''
        return (-np.exp(0.5*self.param_h) + np.exp(-self.param_h*(w-0.5))) / \
              ((1-np.exp(0.5*self.param_h))*(1+np.exp(-self.param_h*(w-0.5))))\
                - self.param_gamma*x

    def _get_system(self, off_nodes=[]):
        '''Build the `dxdt(t, x_array, *args)` right-hand-side function.

        Parameters
        ----------
        off_nodes : iterable of int, optional
            Node **indices** whose derivative should be forced to 0
            (held constant).

        Returns
        -------
        dxdt : callable
            Function ``dxdt(t, x_array, *args) -> ndarray`` computing
            `_dxdt_transform(x, _omega(x))` element-wise, after clipping
            `x_array` to [0, 1]; suitable as the `fun` argument of
            `scipy.integrate.solve_ivp`.
        '''

        def dxdt(t, x_array, *args):
            x_array[x_array < 0] = 0
            x_array[x_array > 1] = 1

            w = self._omega(x_array)
            d = self._dxdt_transform(x_array, w)
            for i in off_nodes:
                d[i] = 0
            return d

        return dxdt


ode_classes = {
    'squad': SquadODE,
    'hillcube': BooleCubeODE,
    'normalisedhillcube': BooleCubeODE,
    'boolecube': BooleCubeODE
}
''' dict : transform to ODE class translation
'''

transforms = set(ode_classes)
''' set : list of accepted ODE transforms
'''
