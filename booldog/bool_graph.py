import os
import numpy as np

from itertools import product
from collections import defaultdict

from pathlib import Path

import warnings

from pyboolnet.trap_spaces import steady_states
from pyboolnet.interaction_graphs import primes2igraph
from pyboolnet.state_transition_graphs import primes2stg, stg2image


from booldog.utils.utils import *


from booldog.simulation_result import BooleanSimulationResult

class BooleanGraph:
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

        logger.info(f"Created Boolean graph with {self.n} nodes")


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
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        intgraph = primes2igraph(self.primes)
        for source, d in intgraph.adjacency(): #nx 2.x
            for target, sub_d in d.items():
                sign = next(iter(sub_d["sign"]))
                if sign == 1:
                    Act[self.index[target], self.index[source]] = 1
                elif sign == -1:
                    Inh[self.index[target], self.index[source]] = 1
                else:
                    print("Warning: Issue with edge: ", target, source)

        return ensure_ndarray(Act), ensure_ndarray(Inh)

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

        else:
            state_array = np.zeros(self.n)
            # these are the fixed points
            for node, state in fixed.items():
                    state_array[self.index[node]] = state

            # generate the free points
            free_variables = list(set(self.nodes) - set(fixed.keys()))
            num_free_variables = len(free_variables)
            print(num_free_variables)
            for free_variable_state in product([0, 1],
                                               repeat=num_free_variables):

                this_state_array = state_array.copy()
                for node, state in zip(free_variables, free_variable_state):
                    this_state_array[self.index[node]] = state
                yield this_state_array

    def boolean_simulation(self,
           initial_values=None,
           plot=True,
           export=False,
           **kwargs):
        '''Compute a Boolean simulation (or state transition)
        from optional initial values.

        Parameters
        ----------
        initial_values : int, str, list, dict, optional
            Initial state, see Notes for format
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
        if initial_values is None:
            stg = primes2stg(self.primes, "synchronous")
        else:
            stg = primes2stg(self.primes, "synchronous", initial_values)

        return BooleanSimulationResult(self, stg)

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
        intgraph = InteractionGraphs.primes2igraph(self.primes)

        if direction == 'out':
            for node, d in intgraph.adjacency(): #nx 2.x
                for other_node, sub_d in d.items():
                    sign = next(iter(sub_d["sign"]))
                    interactions[node][other_node] = sign
        elif direction == 'in':
            for node, d in intgraph.adjacency(): #nx 2.x
                for other_node, sub_d in d.items():
                    sign = next(iter(sub_d["sign"]))
                    interactions[other_node][node] = sign
        else:
            print(f"'Direction' should be 'in' or 'out', not {direction}")

        return interactions



    def __repr__(self):
        #TODO
        return "BooleanGraph()"

    def __len__(self):
        return self.n
