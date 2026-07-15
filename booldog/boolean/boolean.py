'''Module containing functions to manipulate and simulate
Boolean networks. It is not intended to be used directly, but rather as a mixin
for other classes.
'''

from collections import defaultdict
from itertools import product
import re
import os
import logging
import numpy as np

from pyboolnet.interaction_graphs import primes2igraph
from pyboolnet.prime_implicants import is_constant, find_successors, find_predecessors
from pyboolnet.state_space import list_states_in_subspace
from pyboolnet.state_space import state2str
from pyboolnet.state_transition_graphs import primes2stg
from pyboolnet.trap_spaces import compute_steady_states
from pyboolnet.prime_implicants import find_inputs
from pyboolnet.external.bnet2primes import bnet_text2primes

from booldog.simulation_result import BooleanSimulationResult, BooleanStateSpace
from booldog.classes import BoolDogNode
from booldog.utils import ensure_ndarray
from booldog.utils.decorators import validate_node_argument

logger = logging.getLogger(__name__)

class BooleanNetworkMixin():
    '''Class for Boolean network functions.

    This class is not intended to be used directly, but rather as a mixin.
    '''

    def get_primes(self):
        '''Compute the prime implicants of the network from its current node rules.

        Converts the network to bnet format (via :meth:`to_bnet`) and calls
        ``pyboolnet.external.bnet2primes.bnet_text2primes`` on the result.

        Returns
        -------
        primes : dict
            Prime implicants of the Boolean network.

        Raises
        ------
        ValueError
            If the underlying call to bnet2primes fails (returns ``None``),
            e.g. because the bnet text could not be parsed.
        '''

        bnet = self.to_bnet()

        primes = bnet_text2primes(bnet)

        if primes is None:
            raise ValueError("Could not convert bnet to primes.")

        return primes

    @validate_node_argument
    def get_rule(self, node_id):
        '''Get the rule for a given node.

        Parameters
        ----------
        node_id : str or BoolDogNode
            Node identifier or `BoolDogNode` object. Validated and normalised
            to its identifier by the `validate_node_argument` decorator.

        Returns
        -------
        str
            The rule for the given node.

        '''
        if node_id not in self.node_ids:
            raise ValueError(f"{node_id} is not a node in the network.")
        return self.nodes[node_id].rule

    def primes_to_matrices(self) -> tuple[np.ndarray, np.ndarray]:
        '''Represent the Boolean network as an "activation" and an
        "inhibition" matrix, derived from the sign of each edge of the
        network's interaction graph (`pyboolnet.interaction_graphs.
        primes2igraph`).

        Returns
        ----------
        act : np.array
            n * n matrix with entry m_{i, j} = 1 iff node_j activates node_i

        inh : np.array
            n * n matrix with entry m_{i, j} = 1 iff node_j inhibits node_i

        Notes
        ----------
        This only makes sense for networks whose regulatory logic can be
        interpreted as a combination of independent activators/inhibitors
        (i.e. threshold/SQUAD-style networks); it is used by
        `booldog.continuous.ode_factory.SquadODE` to build its ODE system.
        The sign of each edge is taken from the interaction graph built by
        `pyboolnet.interaction_graphs.primes2igraph`; if an edge's sign set
        contains neither 1 nor -1 a warning is logged and the edge is
        otherwise ignored.

        TODO
        ----
        these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to
        use e.g. dok_matrix
        '''
        act_matrix = np.zeros((self.n, self.n))  # activators
        inh_matrix = np.zeros((self.n, self.n))  # inhibitors

        intgraph = primes2igraph(self.primes)
        for source, d in intgraph.adjacency():  #nx 2.x
            for target, sub_d in d.items():
                sign = next(iter(sub_d["sign"]))
                if sign == 1:
                    act_matrix[self.index[target], self.index[source]] = 1
                elif sign == -1:
                    inh_matrix[self.index[target], self.index[source]] = 1
                else:
                    logger.warning("Warning: Issue with edge: %s %s", target,
                                   source)

        return ensure_ndarray(act_matrix), ensure_ndarray(inh_matrix)

    def generate_states(self, fixed=None):
        '''Generate all possible states of the graph.

        Parameters
        ----------
        fixed : dict, optional
            A dictionary of {node: state} for nodes to be kept fixed, with
            node in `self.nodes` and state in {0, 1}. All other nodes are
            free and every combination of their values is generated. If not
            given, every node is free (i.e. all 2**n states are generated).

        Returns
        -------
        iterator
            If `fixed` is None, the cartesian product
            `itertools.product([0, 1], repeat=self.n)` of length-n 0/1
            tuples (all 2**n states, unmapped to node identifiers).
            If `fixed` is given, a generator is returned instead (see
            Yields).

        Yields
        ----------
        state : np.array
            length n array with a state of the graph, indexed by
            `self.index`; only produced when `fixed` is given.

        '''
        if fixed is None:
            # all states = cartesion product
            return product([0, 1], repeat=self.n)

        state_array = np.zeros(self.n)
        # these are the fixed points
        for node, state in fixed.items():
            state_array[self.index[node]] = state

        # generate the free points
        free_variables = list(set(self.nodes) - set(fixed.keys()))
        num_free_variables = len(free_variables)
        logger.debug("Free variables: %i", num_free_variables)
        for free_variable_state in product([0, 1], repeat=num_free_variables):

            this_state_array = state_array.copy()
            for node, state in zip(free_variables, free_variable_state):
                this_state_array[self.index[node]] = state
            yield this_state_array

    def inactivate_state(self):
        '''Build a single-state `BooleanStateSpace` with all nodes inactive
        (i.e. the all-0 state).

        Returns
        -------
        BooleanStateSpace
        '''
        return BooleanStateSpace(self, ["0" * self.n])

    def activate_state(self):
        '''Build a single-state `BooleanStateSpace` with all nodes active
        (i.e. the all-1 state).

        Returns
        -------
        BooleanStateSpace
        '''
        return BooleanStateSpace(self, ["1" * self.n])

    def boolean_simulation(self, initial_states=None):
        '''Compute a Boolean simulation (or state transition)
        from optional initial values.

        Parameters
        ----------
        initial_states : str, list, dict, callable, or BooleanStateSpace, optional
            Initial states, see Notes for format. If not given, every state
            is treated as initial (the full state transition graph is
            built).

        Returns
        -------
        BooleanSimulationResult
            The (synchronous) state transition graph, together with the
            standardised initial states.

        Notes
        ----------
        This is a wrapper for pyboolnet.state_transition_graphs.primes2stg
        with the update scheme fixed to ``"synchronous"``, and therefore
        takes the same argument format for initial states.

        **From pyboolnet documentation:**

        .. code-block:: text

            Either a list of states in dict or str format::

                init = ["000", "111"]
                init = ["000", {"v1":1,"v2":1,"v3":1}]

            or as a function that is called on every state and must return
            either True or False to indicate whether the state ought to be initial::

                init = lambda x: x["v1"]>=x["v2"]

            or by a subspace in which case all the states contained in it are
            initial::

                init = "--1"
                init = {"v3":1}
        '''
        if initial_states is None:
            stg = primes2stg(self.primes, "synchronous")
        else:

            if isinstance(initial_states, BooleanStateSpace):
                initial_states = initial_states.state_space

            stg = primes2stg(self.primes, "synchronous", initial_states)

        return BooleanSimulationResult(
            self, stg, self.standard_states_format(initial_states))

    def standard_states_format(self, states):
        '''Convert a set of states, given in any of the formats accepted by
        `boolean_simulation`/`pyboolnet.state_transition_graphs.primes2stg`,
        to a plain list of state strings.

        Parameters
        ----------
        states : None, str, dict, list, or callable
            - None: returns None (no standardisation is performed).
            - a callable taking a state in dict format (`{node_id: 0/1}`)
              and returning a bool: every state of the full state space for
              which it returns True is included.
            - a str or dict subspace (e.g. ``"--1"`` or ``{"v3": 1}``): all
              states contained in that subspace are included (via
              `pyboolnet.state_space.list_states_in_subspace`).
            - otherwise (e.g. a list of str/dict states): each element is
              converted to its string representation.

        Returns
        -------
        list of str, or None
            The states in string representation, or None if `states` is
            None.

        Notes
        -----
        Mirrors the state-format handling of pyboolnet's ``primes2stg``, see
        https://github.com/hklarner/pyboolnet/blob/529860bc1185277fb2b5e0f3b36c9ba6c7b9fe2f/pyboolnet/state_transition_graphs.py#L124-L135
        '''

        if states is None:
            return None

        space = len(self.nodes) * [[0, 1]]

        # standardise the initial states
        if hasattr(states, "__call__"):
            fringe = [
                dict(zip(self.nodes, values)) for values in product(*space)
            ]
            fringe = [state2str(x) for x in fringe if states(x)]

        elif type(states) in [str, dict]:
            fringe = list_states_in_subspace(primes=self.primes,
                                             subspace=states)

        else:
            fringe = [state2str(x) for x in states]

        return fringe

    def steady_states(self):
        '''Compute all steady states (fixed points) of the Boolean network.

        Returns
        -------
        BooleanStateSpace
            The steady states of the network.

        Notes
        -----
        This is a wrapper for `pyboolnet.trap_spaces.compute_steady_states`,
        which uses the Potassco ASP solver to compute steady states as a
        special case of trap spaces (using its defaults, at most 1000
        steady states are returned).
        '''

        all_steady_states = compute_steady_states(self.primes)
        return BooleanStateSpace(self, all_steady_states)

    @validate_node_argument
    def get_parents(self, node_id):
        '''Fetch regulators/inputs to a node.

        Parameters
        ----------
        node_id : str or BoolDogNode
            Node identifier or `BoolDogNode` object. Validated and
            normalised to its identifier by the `validate_node_argument`
            decorator.

        Returns
        ----------
        parents : list
            Sorted list of identifiers of the nodes that regulate `node_id`
            (i.e. appear in its update rule).

        Notes
        -----
        This is a wrapper for `pyboolnet.prime_implicants.find_predecessors`.
        '''
        return find_predecessors(self.primes, [node_id])

        # set([key for d in self.primes[node][0] for key in d.keys()] +\
                   # [key for d in self.primes[node][1] for key in d.keys()])

    def get_interactions(self, direction='out'):
        '''Build a nested dict of pairwise regulatory signs from the
        network's interaction graph.

        Parameters
        ----------
        direction : {'out', 'in'}, optional
            direction out: interactions[a][b] describes the edge a --> b
            (`a` regulates `b`).
            direction in: interactions[a][b] describes the edge a <-- b
            (`b` regulates `a`).
            Any other value logs a warning and returns an empty dict.

        Returns
        -------
        interactions : collections.defaultdict of dict
            interactions[a][b] = sign, where sign is 1 (activation) or -1
            (inhibition) — one of the signs from the interaction-graph edge
            (see `pyboolnet.interaction_graphs.primes2igraph`); an edge with
            both signs recorded (dual interaction) is not distinguished
            from a single-signed edge.

        Notes
        -----
        Wraps `pyboolnet.interaction_graphs.primes2igraph`. See also
        `primes_to_matrices`, which builds a similar (but dense-matrix,
        binary) representation.
        '''
        interactions = defaultdict(dict)
        intgraph = primes2igraph(self.primes)

        # TODO this is wrong?
        if direction == 'out':
            for node, d in intgraph.adjacency():  #nx 2.x
                for other_node, sub_d in d.items():
                    sign = next(iter(sub_d["sign"]))
                    interactions[node][other_node] = sign
        elif direction == 'in':
            for node, d in intgraph.adjacency():  #nx 2.x
                for other_node, sub_d in d.items():
                    sign = next(iter(sub_d["sign"]))
                    interactions[other_node][node] = sign
        else:
            logger.warning('"direction" should be "in" or "out", not "%s"',
                           direction)

        return interactions

    @validate_node_argument
    def is_constant(self, node_id):
        '''
        Whether node is a constant (i.e. its rule is a fixed 0 or 1,
        independent of any other node) in the network.

        Parameters
        ----------
        node_id : str or BoolDogNode
            Node identifier or `BoolDogNode` object. Validated and
            normalised to its identifier by the `validate_node_argument`
            decorator.

        Returns
        ----------
        constant : bool
            Whether `node_id` is a constant.

        Notes
        -----
        This is a wrapper for `pyboolnet.prime_implicants.is_constant`.

        '''

        return is_constant(self.primes, node_id)

    def list_network_inputs(self):
        '''
        List all input nodes of the network, i.e. nodes whose own value is
        their only regulator (`node, node` as a bnet rule) so that their
        value never changes under simulation.

        Returns
        ----------
        inputs : list
            Sorted list of input node identifiers.

        Notes
        -----
        This is a wrapper for `pyboolnet.prime_implicants.find_inputs`.

        '''

        inputs = find_inputs(self.primes)
        return inputs
