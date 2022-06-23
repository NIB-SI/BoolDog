import numpy as np
from itertools import product

import logging
logger = logging.getLogger(__name__)

from booldog.utils.utils import *

##############################
#      GENERAL FUNCTIONS     #
##############################



##############################
#        CHILD CLASS         #
##############################

def ODE_factory(graph, transform, **kwargs):
    '''Create a :py:class:`booldog.ode.ODE` from :py:class:`RegulatoryNetwork`.

    Parameters
    ----------
    graph : :py:class:`RegulatoryNetwork` or :py:class:`BooleanGraph`
        Input graph.

    transform : str
        One of accepted transforms. See :py:data:`booldog.ode.transforms` or :ref:`Notes <tagnotesode>` for options.

    Other Parameters
    ----------------
    **kwargs
        Additional arguments and keyword arguments passed to specific parent
        class initializer.

    Returns
    -------
    ode : :py:class:`booldog.ode.ODE`
        A ODE system


    .. _tagnotesode:

    Notes
    -----
    For specific transforms, see the relevant parent class for keyword
    arguments (`**kwargs`).

    The parent class per transform is defined in
    :py:data:`booldog.ode.ode_parent_classes`.


    If the parameter is an int or float, the value is assigned for
    all variables. Otherwise the parameter argument should be a dict
    with keys as node names and values for their initial state. In
    this case, if the initial state is not defined for all nodes, a
    `default` key with the default value should also be present in the
    dict.

    Here follows a summary of the transform-specific keyword arguments

        'squad'

            - :py:class:`booldog.ode.SquadODE`
            - gamma : self-decay
            - h :  sigmoid gain

        'boolecube'

            - :py:class:`booldog.ode.BoolCubeODE`
            - tau : life-time of species

        'hillcube'

            - :py:class:`booldog.ode.BoolCubeODE`
            - tau : life-time of species
            - n : Hill coefficient
            - k : Hill dissociation constant

        'normalisedhillcube'

            - :py:class:`booldog.ode.BoolCubeODE`
            - tau : life-time of species
            - n : Hill coefficient
            - k : Hill dissociation constant
    '''
    transform = transform.lower()
    if transform == 'placeholder':
        ParentClass = BoolCubeODE
    elif not transform in ode_parent_classes.keys():
        raise ValueError(f"transform' argument must be one of"\
                         f"{list(ode_parent_classes.keys())}")
    else:
        ParentClass = ode_parent_classes[transform]

    logger.info(f"Creating {ParentClass} ODE object")

    class ODE(ParentClass):
        '''Generic ODE class produced by factory.

        Parent class is a variable, and defined by `transform` argument. '''

        def __init__(self, graph, transform, **kwargs):
            '''Initialise ODE

            Parameters
            ----------
            graph : :py:class:`booldog.BooleanGraph`
            transform : str

            Other Parameters
            ----------
            **kwargs
                In the case `graph` is not a BooleanGraph instance, additional
                keyword arguments passed to :py:class:`booldog.BooleanGraph`.
            '''
            if transform == 'placeholder':
                return

            # if not isinstance(graph, RegulatoryNetwork):
            #     raise TypeError(f"'graph' argument must be a RegulatoryNetwork object."\
            #                     f"not {type(graph)}. ")

            parent_class_name = self.__class__.__mro__[1].__name__
            logger.info(f"Initialising ODE system ... with class {parent_class_name}. ")

            self.n = len(graph)
            self.boolean_graph = graph
            self.transform = transform

            super().__init__(transform, **kwargs)

            logger.info("Done. ")


        def event_function(self, t, x, event_t, *args):
            '''Event function for `events` of `scipy.integrate.solve_ivp`.

            Parameters
            ----------
            t : float
                Time-point of simulation
            x : narray

            event_t : float
                Time-point of event

            *args
                ignored

            Attributes
            ----------
            terminal : True

            '''
            return t - event_t
        event_function.terminal = True

        def update(self, off_nodes=None):
            ''' Resets dxdt

            Shortcut to parent class's `_get_system` method.

            Parameters
            ----------
            off_nodes : list of int, optional
                List of node **indices** to set derivative to 0,
                i.e. these nodes will remain constant in simulation.
            '''
            self.dxdt =  self._get_system(off_nodes=off_nodes)

    # for sphinx documentation
    ODE_factory.ex_class = ODE
    ODE_factory.ex_class.__bases__ = tuple(set(ode_parent_classes.values()))

    return ODE(graph, transform, **kwargs)

##############################
#       PARENT CLASSES       #
##############################

# https://github.com/krumsieklab/Odefy/blob/11d048d550a8f64250ba01f76f5a83048c8be6cf/Odefy-1.20/code/models/CreateCubeCalls.m
class BoolCubeODE():
    '''An ODE parent class.

    Use of multivariate polynomial interpolation for the transformation of a Boolean graph to a system of ODEs.

    Attributes
    ----------
    dxdt : function

    param_tau : arraylike
        life-time of species

    param_n : arraylike
        Hill coefficient

    param_k : arraylike
        Hill dissociation constant

    param_dict : dict
        track parameters
     '''

    def __init__(self, transform, tau=1, n=3, k=0.5, **kwargs):
        ''' Initialise BoolCube ODE system.

        Parameters
        ----------
        transform : str

        tau : int, float, or dict, optional
            life-time of species

        n : int, float, or dict, optional
            Hill coefficient

        k : int, float, or dict, optional
            Hill dissociation constant

        verbose : bool

        Notes
        ----------
        Only used as Parent class.

        tau_i = zero --> dx_i/dt = 0

        References
        ----------
        [1] Wittmann, D. M., Krumsiek, J., Saez-Rodriguez, J., Lauffenburger,
        D., A., Klamt, S., & Theis, F. J. (2009). Transforming Boolean models
        to continuous models: Methodology and application to T-cell receptor
        signaling. BMC Systems Biology, 3(1), 98.
        https://doi.org/10.1186/1752-0509-3-98

        '''

        self.param_n = parameter_to_array(n, self.boolean_graph.index)
        self.param_k = parameter_to_array(k, self.boolean_graph.index)

        self.param_tau = parameter_to_array(tau, self.boolean_graph.index)

        self.param_dict = {"n":self.param_n,
                           "k":self.param_k,
                           "tau":self.param_tau}

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

        # returns an array function
        self.dxdt =  self._get_system()

    def hill(self, x_array):
        return x_array**self.param_n / \
               (x_array**self.param_n + self.param_k**self.param_n)

    def normalised_hill(self, x_array):
        return self.hill(x_array) / self.hill(1)

    def identity(self, x_array):
        return x_array

    def _get_system(self, off_nodes=None):

        if off_nodes is None:
            off_nodes = set()
        else:
            off_nodes = set(off_nodes)

        off_nodes.update(np.where(self.param_tau ==0)[0])

        def dxdt(t, x_array, *args):
            x_array[x_array < 0] = 0
            x_array[x_array > 1] = 1

            b = self.B1(self.transform_function(x_array))
            d = 1/self.param_tau * ( b  - x_array)
            for i in off_nodes:
                d[i] = 0
            return d

        return dxdt

    def homologue_b1(self):
        ''' Create function to calculate the multivariate polynomial
        interpolation of Boolean functions

        Returns
        ----------
        B1 : function

        '''
        # spaces = set()
        # sums = []
        # all_B1s = ['0']*self.boolean_graph.n
        # for node in self.boolean_graph.nodes: # iterate over all nodes
        #     for prime_dict in self.boolean_graph.primes[node][1]:
        #         for x_bool in self.boolean_graph.generate_states(
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
        #         all_B1s[self.boolean_graph.index[node]]  = B1

        all_B1s = ['0']*self.boolean_graph.n
        for node in self.boolean_graph.nodes: # iterate over all nodes
            parents = self.boolean_graph.get_parents(node)

            states = []
            spaces = set()
            for prime_dict in self.boolean_graph.primes[node][1]:

                fixed_parents = prime_dict.keys()
                free_parents = parents - set(fixed_parents)

                if len(free_parents) > 0:
                    for x in product([0, 1], repeat=len(free_parents)):
                        this_state = prime_dict.copy()
                        for parent, parent_state in zip(free_parents, x):
                            this_state[parent] = parent_state
                        str_rep = "".join([str(this_state[node]) if node in this_state else "-" for node in self.boolean_graph.nodes])
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
                str_subterms  = []
                for this_node, this_node_state in state.items():
                    if this_node_state == 0:
                        subterms.append(f'(1-x[{self.boolean_graph.index[this_node]}])')
                        str_subterms.append(f'(1-x[{this_node}])')
                    else:
                        subterms.append(f'x[{self.boolean_graph.index[this_node]}]')
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
                    print(f"The rule for node {node} is too long (depends on {len(states)} states. ")

                all_B1s[self.boolean_graph.index[node]]  = B1

            logger.debug('%s %i', node, len(states))


        return eval('lambda x:' + 'np.array([' + ','.join(all_B1s) + '])')
        #return lambda x: x

    def write_c_code():
        pass


# class SquadODE():
#     '''An ODE parent class.

#     Use of SQUAD for the transformation of a Boolean graph to a system of ODEs.

#     Attributes
#     ----------
#     dxdt : function

#     param_gamma : arraylike
#         decay rate

#     param_h : arraylike
#         sigmoidal gain

#     Act : arraylike
#         activator matrix

#     Inh : arraylike
#         inhibitor matrix

#     '''

#     def __init__(self, transform, gamma=1, h=10, **kwargs):
#         '''Transform a Activation/Inhibition Boolean network
#         into an ODE system via SQUAD transform.

#         Parameters
#         ----------

#         transform : str

#         gamma : int, float, or dict, optional
#             decay rate

#         h : int, float, or dict, optional
#             sigmoidal gain

#         Notes
#         ----------
#         Only used as Parent class.

#         References
#         ----------
#         [1] Di Cara, A., Garg, A., De Micheli, G., Xenarios, I., & Mendoza, L.
#         (2007). Dynamic simulation of regulatory networks using SQUAD.
#         BMC Bioinformatics, 8(1), 1–10. https://doi.org/10.1186/1471-2105-8-462
#         '''

#         self.param_gamma = parameter_to_array(gamma, self.boolean_graph.index)
#         self.param_h = parameter_to_array(h, self.boolean_graph.index)

#         print(self.param_gamma)
#         print(self.param_h)
#         # matrices
#         self.Act, self.Inh = self.boolean_graph.primes_to_matrices()

#         # needed for computations
#         col_ones = np.ones((self.n))
#         self._A1 = self.Act.dot(col_ones)
#         self._a1 = (1+self._A1)/self._A1
#         self._B1 = self.Inh.dot(col_ones)
#         self._b1 = (1+self._B1)/self._B1

#         self.dxdt = self._get_system()

#     def _omega(self, x):
#         '''
#         Equation (2) of Di Cara et al (2007)
#         http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-462
#         Based on Andre Blejec R code.
#         '''

#         x = ensure_ndarray(x)

#         col_ones = np.ones((self.n))

#         Ax = self.Act.dot(x)
#         a = ensure_ndarray( self._a1 * Ax/(1+Ax) )
#         a[~np.isfinite(a)] = 1

#         Bx = self.Inh.dot(x)
#         b = ensure_ndarray( self._b1 * Bx/(1+Bx) )
#         b[~np.isfinite(b)] = 0

#         o = ensure_ndarray(a * (1-b))
#         o[np.where(self._A1 + self._B1 ==0)] = 0

#         return o

#     def _dxdt_transform(self, x, w):
#         ''' Equation (2) of Di Cara et al (2007) '''
#         return (-np.exp(0.5*self.param_h) + np.exp(-self.param_h*(w-0.5))) / \
#               ((1-np.exp(0.5*self.param_h))*(1+np.exp(-self.param_h*(w-0.5))))\
#                 - self.param_gamma*x

#     def _get_system(self, off_nodes=[]):

#         def dxdt(t, x_array, *args):
#             x_array[x_array < 0] = 0
#             x_array[x_array > 1] = 1

#             w = self._omega(x_array)
#             d = self._dxdt_transform(x_array, w)
#             for i in off_nodes:
#                 d[i] = 0
#             return d

#         return dxdt


# class ShaoODE():
#     '''
#     Source: From Boolean Network Model to Continuous Model Helps in Design of
#     Functional Circuits
#     '''
#     def __init__(self, graph, **kwargs):
#         super().__init__(graph)
#         self.dxdt = self._squad(graph, gamma, h)





# class RacipeODE():

#     def __init__(self, graph, transform, gamma=1, h=10, **kwargs):
#         '''
#         Transform a Activation/Inhibition Boolean network
#         into an ODE system via the RACIPE transform.

#         If ODE parameters are passed as an int or float, the value is
#         assigned for all variables
#         other wise a dict with a 'default' key and value, and key value
#         pairs for all others.

#         Parameters
#         ----------

#         gamma : float, dict
#             decay rate

#         h : float, dict
#             gain(?)
#         '''
#         super().__init__(graph)

#         self.dxdt = self._squad(graph, gamma, h)
#         print("done. ")




ode_parent_classes = {
    # 'squad':SquadODE,
    # 'shao':ShaoODE,
    'hillcube':BoolCubeODE,
    'normalisedhillcube':BoolCubeODE,
    'boolecube':BoolCubeODE
}
''' dict : transform to parent class translation
'''

transforms = set([
    'boolecube',
    'hillcube',
    'normalisedhillcube',
    # 'shao',
    # 'squad'
])
''' set : list of accepted ODE transforms
'''
