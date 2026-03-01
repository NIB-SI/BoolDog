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
        ''' Get prime implicants
        Parameters
        ----------

        Returns
        -------
        primes : dict
            Prime implicants of the Boolean network.
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
        node : str or BoolDogNode
            Node identifier or BoolDogNode object.

        Returns
        -------
        str
            The rule for the given node.

        '''
        if node_id not in self.node_ids:
            raise ValueError(f"{node_id} is not a node in the network.")
        return self.nodes[node_id].rule

    def primes_to_matrices(self) -> tuple[np.ndarray, np.ndarray]:
        '''Represent Boolean network a "activation" and "inhibition" matrices.

        Returns
        ----------
        act : np.array
            n * n matrix with entry m_{i, j} = 1 iff node_j activates node_i

        inh : np.array
            n * n matrix with entry m_{i, j} = 1 iff node_j inhibits node_i

        Notes
        ----------
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Used in SquadODE
        Create logic matrices

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
        fixed : dict
            A dictionary of {node:state} to be kept fixed, with
            node in graph.nodes, and state in {0, 1}.

        Yields
        ----------
        state : np.array
            length n array with a state of the graph.

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
        '''A state space with all nodes inactive'''
        return BooleanStateSpace(self, ["0" * self.n])

    def activate_state(self):
        '''A state space with all nodes active'''
        return BooleanStateSpace(self, ["1" * self.n])

    def boolean_simulation(self, initial_states=None):
        '''Compute a Boolean simulation (or state transition)
        from optional initial values.

        Parameters
        ----------
        initial_states : int, str, list, dict, or BooleanStateSpace, optional
            Initial states, see Notes for format

        Notes
        ----------
        This is a wrapper for pyboolnet.state_transition_graphs.primes2stg,
        and therefore takes the same argument format for initial states.

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
        '''

        Notes
        -----

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
        '''All steady states of the Boolean graph.
        '''

        all_steady_states = compute_steady_states(self.primes)
        return BooleanStateSpace(self, all_steady_states)

    @validate_node_argument
    def get_parents(self, node_id):
        '''Fetch regulators/inputs to a node

        Parameters
        ----------
        node : str
            Node name

        Returns
        ----------
        parents : set
            Set of parent nodes

        '''
        return find_predecessors(self.primes, [node_id])

        # set([key for d in self.primes[node][0] for key in d.keys()] +\
                   # [key for d in self.primes[node][1] for key in d.keys()])

    def get_interactions(self, direction='out'):
        '''
        direction out: interactions[a][b]  a --> b
        direction in:  interactions[a][b]  a <-- b
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
        Whether node is a constant (has no input) in the network.

        Parameters
        ----------
        node_id : str
            Node identifier

        Returns
        ----------
        constant : bool
            Whether `node` is a constant

        Notes
        -----
        This is a wrapper for pyboolnet.prime_implicants.is_constant.

        '''

        return is_constant(self.primes, node_id)

    def list_network_inputs(self):
        '''
        List all nodes that have no regulators/inputs in the network.

        Returns
        ----------
        inputs : list
            List of input node names

        Notes
        -----

        '''

        inputs = find_inputs(self.primes)
        return inputs
