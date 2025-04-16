from collections import defaultdict
from itertools import product

import logging
import numpy as np
from enum import Enum

from pyboolnet.interaction_graphs import primes2igraph
from pyboolnet.prime_implicants import is_constant, find_successors, create_variables, remove_variables
from pyboolnet.state_space import list_states_in_subspace
from pyboolnet.state_space import state2str
from pyboolnet.state_transition_graphs import primes2stg
from pyboolnet.trap_spaces import steady_states

from booldog.simulation_result import BooleanSimulationResult
from booldog.utils import ensure_ndarray, ExtendedEnum

logger = logging.getLogger(__name__)


class ModificationTypes(ExtendedEnum):
    ADD = "add_node"
    REMOVE = "remove_node"
    UPDATE = "update"

class BooleanNetworkModification():

    def __init__(self, modification_type, node, rule=None):
        self.type = modification_type
        self.node = node
        self.rule = rule


class BooleanNetwork():
    '''A class to represent a Boolean graph.

    Attributes
    ----------
    n : int
        The number of nodes/variables in the graph

   primes : dict
        Prime implicants of the Boolean graph. See
        `PyBoolNet:prime implicants
        <https://pyboolnet.readthedocs.io/en/latest/Manual.html#prime-implicants>`_
        for more information.

    nodes : tuple of str
        Lists node names in the graph

    index : dict
        Dictionary of node name to integer index for indexing arrays
    '''

    def __init__(self, primes, nodes, index, n):
        '''Initialise a Boolean graph.

        Parameters
        ----------
        graph : str or dict
            A file path to the graph or a dictionary with the graph.

        data_format : {'primes', 'interactions', 'bnet', 'graphml'}
            String specifying data format.

        kwargs
            #TODO Additional keyword arguments for the importer function.

        '''
        self.primes, self.nodes, self.index, self.n = primes, nodes, index, n

        self._modifications = []

    def primes_to_matrices(self):
        '''Reduce graph to Activation and Inhibition Matrices.

        Returns
        ----------
        Act : np.array
            n * n matrix with entry m_{i, j} = 1 iff node_j activates node_i

        Inh : np.array
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

    def boolean_simulation(
        self,
        initial_states=None
    ):
        '''Compute a Boolean simulation (or state transition)
        from optional initial values.

        Parameters
        ----------
        initial_states : int, str, list, dict, optional
            Initial states, see Notes for format
        plot : bool, optional
            Whether to plot the simulation results
        export : bool or path object or string, optional
            False, or a file path to save simulation results.
            Exports valu

        Other Parameters
        ----------------
        **kwargs
            If `plot=True` , additional keyword arguments are
            passed to :py:func:`plot_boolean_simulation`.

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

        all_steady_states = steady_states(self.primes)
        return all_steady_states

    def get_parents(self, node):
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

        return set([key for d in self.primes[node][0] for key in d.keys()] +\
                   [key for d in self.primes[node][1] for key in d.keys()])

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

    def is_constant(self, node):
        '''
        Whether node is a constant (has no input) in the network.

        Parameters
        ----------
        node : str
            Name of node

        Returns
        ----------
        constant : bool
            Whether `node` is a constant

        Notes
        -----
        This is a wrapper for pyboolnet.prime_implicants.is_constant.

        '''

        return is_constant(self.primes, node)

    #############################
    #    EDITING THE NETWORK    #
    #############################

    def _reset_from_primes(self):
        '''Helper function to reset attributes after primes have been updated. '''
        self.nodes = tuple(sorted(
            self.primes.keys()))  # tuple (i.e. immutable)
        self.index = {i: node for node, i in enumerate(self.nodes)}
        self.n = len(self.nodes)


    def _track_network_modifications(self, modification_type, node, rule=None):
        '''
        Keep track of network modifications
        '''
        self._modifications.append(BooleanNetworkModification(modification_type, node, rule=rule))

    def modify(self, modifications):

        for modification in modifications:
            match modification.modification_type:
                case ModificationTypes.ADD:
                    self.add_node(modification.node, modification.rule)
                case ModificationTypes.REMOVE:
                    self.remove_nodes(modification.node)
                case ModificationTypes.UPDATE:
                    self.update_rule(modification.node, modification.rule)
                case _:
                    raise ValueError(f"Modification type has to be one of {', '.join(ModificationTypes.values())}, not {modification.modification_type:}")

    def remove_node(self, node):
        '''
        Removes `node` from the network.

        Parameters
        ----------
        node : str
            Names of node to remove

        Returns
        ----------
        None

        Notes
        -----
        A node cannot be removed if other nodes depend on it (i.e. it occurs in
        their update logic). To remove such a node, either also remove all of
        its dependants, or first update the logic rule of its dependants to
        remove dependency.

        This is a wrapper for pyboolnet.prime_implicants.remove_variables.
        '''

        self.remove_nodes(node)

    def remove_nodes(self, nodes):
        '''
        Removes all nodes in `nodes` from the network.

        Parameters
        ----------
        nodes : list
            List of names of nodes to remove

        Returns
        ----------
        None

        Notes
        -----
        A node cannot be removed if other nodes depend on it (i.e. it occurs in
        their update logic). To remove such a node, either also remove all of
        its dependants, or first update the logic rule of its dependants to
        remove dependency.

        This is a wrapper for pyboolnet.prime_implicants.remove_variables.
        '''
        if isinstance(nodes, str):
            nodes = [nodes]

        hit = {
            node: [
                x for x in find_successors(primes=self.primes, sources=[node])
                if x not in nodes
            ]
            for node in nodes
        }
        s = ''
        for node in hit:
            if hit[node]:
                s += f'Cannot remove a node that has dependents. To remove '\
                f'"{node}", you need to remove its dependants as well, or '\
                f'remove their dependency on "{node}" by updating their rules. '\
                f'Dependent(s) for "{node}" are: {", ".join(hit[node])}.\n'
                # TODO example
        if s:
            raise ValueError(s)

        remove_variables(self.primes, nodes)
        self._reset_from_primes()

        self._track_network_modifications(ModificationTypes.REMOVE, nodes)

    def add_node(self, node, rule):
        '''
        Add a new node to Network.

        Parameters
        ----------
        node : str
            Node name
        rule : str
            Rule to define update of `node`, in bnet form. All nodes in
            `rule` need to be defined in Network.

        Returns
        ----------
        None

        This is a wrapper for pyboolnet.prime_implicants.create_variables.
        '''
        if node in self.nodes:
            raise ValueError("f{node} is already present in Network.")

        create_variables(self.primes, {node: rule})
        self._reset_from_primes()

        self._track_network_modifications(ModificationTypes.ADD, node, rule=rule)

    def update_rule(self, node, rule):
        '''
        Update (modify) or add the logic rule defining the update of `node`.
        If `node` does not yet exist, it will be added, if it does exist, its
        update logic will be overwritten.

        Parameters
        ----------
        node : str
            Node name
        rule : str
            New rule to define update of `node`, in bnet form. All nodes in
            `rule` need to be defined in Network.

        Returns
        ----------
        None

        This is a wrapper for pyboolnet.prime_implicants.create_variables.
        '''
        create_variables(self.primes, {node: rule})
        self._reset_from_primes()

        self._track_network_modifications(ModificationTypes.UPDATE, node, rule=rule)

    #############################
    #       REPRESENTATION      #
    #############################

    def __repr__(self):
        return f"{self.__class__} with {self.n} nodes"

    def __len__(self):
        return self.n


        '''
        xxx

        Parameters
        ----------
        x : x
            x

        Returns
        ----------
        x : x
            x

        '''