import os
import numpy as np

from itertools import product
from collections import defaultdict

from pathlib import Path
import importlib

from PyBoolNet import AspSolver
from PyBoolNet import InteractionGraphs
from PyBoolNet import StateTransitionGraphs

from .utils import io
from .utils.utils import *


class BooleanGraph:
    '''
    A class to represent a Boolean graph.

    Attributes
    ----------

    n : int
        The number of nodes/variables in the graph
    primes

    nodes

    index
    '''

    def __init__(self, graph, data_format='bnet', **kwargs):
        '''Initialise a Boolean graph.

        Parameters
        ----------
        graph : str or dict
            A file path to the graph or a dictionary with the graph.
        data_format : str
            String specifying data format.

                * primes
                * interactions
                * bnet
                * graphml

        kwargs #TODO
            Additional keyword arguments for the importer function.

        '''
        if isinstance(graph, BooleanGraph):
            self.primes, self.nodes, self.index, self.n = \
                    graph.primes, graph.nodes, graph.index, graph.n

        else:
            self.primes, self.nodes, self.index, self.n =\
                                    io.read_boolean_graph(graph, data_format, **kwargs)

        print(f"Imported Boolean graph with {self.n} nodes")


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

        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to
        use e.g. dok_matrix
        '''
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        intgraph = InteractionGraphs.primes2igraph(self.primes)
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

    def generate_states(self, fixed={}):
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
        if len(fixed) == 0:
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

    def plot_state_transitions(self, fig,
                               initial_values=None,
                               new_style=True):
        '''Plot the state transition graph, from optional initial values.

        Parameters
        ----------
        fig : str
            File name of generated figure

        initial_values : int, str, list, dict
            Initial state, see Notes for format

        new_style : bool
            Whether to use booldog style, or PyBoolNet style (default) to
            plot the state transition graph. Requires pygraphviz (If not
            installed, will default to PyBoolNet version.)

        Notes
        ----------
        This is a wrapper for PyBoolNet.StateTransitionGraphs.primes2stg,
        and therefore takes the same argument format for initial states.

        From PyBoolNet documentation:
        > Either a list of states in dict or str format

            init = ["000", "111"]
            init = ["000", {"v1":1,"v2":1,"v3":1}]

        > or as a function that is called on every state and must return
        > either True or False to indicate whether the state ought to be initial:

            init = lambda x: x["v1"]>=x["v2"]

        > or by a subspace in which case all the states contained in it are
        > initial:

            init = "--1"
            init = {"v3":1}

        #TODO
        optional list of nodes to plot

        '''
        if initial_values is None:
            stg = StateTransitionGraphs.primes2stg(self.primes,
                                                   "synchronous")
        else:
            stg = StateTransitionGraphs.primes2stg(self.primes,
                                                   "synchronous",
                                                   initial_values)

        if new_style and (importlib.util.find_spec('pygraphviz') is not None):
            import networkx as nx

            rev_index = {}
            for node in self.nodes:
                rev_index[self.index[node]] = node


            stg.graph['node']['shape'] = 'plaintext'
            colors = {"1":"#80ff8a", "0":"#ff9580"}

            for n in stg:
                label = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" ><TR>'
                for i, x in enumerate(n):
                    label += f'<TD BGCOLOR="{colors[x]}">{rev_index[i][:5]} <BR/> {x} </TD>'
                label += """</TR></TABLE>>"""
                stg.nodes[n]["label"]=label

            A = nx.drawing.nx_agraph.to_agraph(stg)
            A.layout('dot')
            A.draw(fig)
            print(f"Saved figure to {fig}. ")

        else:
            if new_style:
                print("pygraphviz not available, defaulting to PyBoolNet")

            stg.graph["node"]["color"] = "cyan"
            stg.graph["node"]["height"] = 0.3
            stg.graph["node"]["width"] = 0.45
            StateTransitionGraphs.stg2image(stg, fig, LayoutEngine="dot")

    def steady_states(self):
        '''Returns all steady states of the Boolean graph.

        '''
        steady_states = AspSolver.steady_states(self.primes)
        return steady_states


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



    def __print__():
        #TODO
        pass

    def __len__(self):
        return self.n
