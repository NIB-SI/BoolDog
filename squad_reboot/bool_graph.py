import os
import numpy as np
from itertools import product

from PyBoolNet import QuineMcCluskey as QMC
from PyBoolNet import AspSolver
from PyBoolNet import InteractionGraphs


from . import io
from .utils import *
from pathlib import Path


def boolean_graph_format_io(graph, in_format, out_format, outfile):
        formats = {'primes':0, 'bnet':0, 'interactions':0,
                            'graphml':0}

        if not ( (in_format in formats) and (out_format in formats)  ):
               raise ValueError("Both 'format' arguments must be in%s"%str(list(formats.keys())))


        if not isinstance(graph, (dict, str, Path)):
               raise TypeError("'graph' argument must be a dictionary or string/file path")

        if isinstance(graph, (str, Path)):
            file_path = Path(graph)
            if not file_path.exists():
               raise FileNotFoundError("'graph' is a path or string, but the file %s does not exist."%graph)

        outfile = Path(outfile)
        file_writable(outfile)



        if (format == 'primes') and isinstance(data, dict):
            # nothing to do
            self.primes = data

        elif (format == 'primes') and (os.path.isfile(data)):
            # primes from a file
            self.primes = io.import_primes(data)

        elif (format == 'interactions') and isinstance(data, dict):
            self._squad_init(data)

        elif (format == 'bnet') and (isinstance(data, str) or os.path.isfile(data)):
            # primes from a file
            self.primes = io.import_bnet(data)

        elif (format == 'graphml') and (os.path.isfile(data)):
            self.primes = io.import_graphml(path, **kwargs)


        else:
            raise NotImplementedError("import of %s from %s is not implemeted. "%(format, str(type(data))))



















class BooleanGraph:
    '''
    A class to represent a Boolean graph.

    Attributes
    ----------

    n : int
        The number of nodes/variables in the graph

    Methods
    ----------

    primes_to_matrices

    interactions_to_matrices

    generate_states

    plot_transitions

    bool_steady_states

    '''

    def __init__(self, data, data_format='bnet', **kwargs):
        '''
        Initialise a Boolean graph.

        Parameters
        ----------

        data : str or dict
            A file path to the graph or a dictionary with the graph.
        data_format : str
            String specifying data format.

                * primes
                * interactions
                * bnet
                * graphml

        kwargs
            Additional keyword arguments for the importer function.

        '''
        format_classes = {'primes':0, 'bnet':0, 'interactions':0,
                            'graphml':0}

        if not data_format in format_classes.keys():
               raise ValueError("'data_format' argument must be one of %s"%str(list(format_classes.keys())))


        if not isinstance(data, (dict, str)):
               raise TypeError("'data' argument must be a dictionary or string")


        if (data_format == 'primes') and isinstance(data, dict):
            # nothing to do
            self.primes = data

        elif (data_format == 'primes') and (os.path.isfile(data)):
            # primes from a file
            self.primes = io.import_primes(data)

        elif (data_format == 'interactions') and isinstance(data, dict):
            self._squad_init(data)

        elif (data_format == 'bnet') and (isinstance(data, str) or os.path.isfile(data)):
            # primes from a file
            self.primes = io.import_bnet(data)

        elif (data_format == 'graphml') and (os.path.isfile(data)):
            self.primes = io.import_graphml(path, **kwargs)


        else:
            raise NotImplementedError("import of %s from %s is not implemeted. "%(data_format, str(type(data))))

        if not (data_format == 'interactions'):
            self.nodes = tuple(sorted(self.primes.keys())) # tuple (i.e. immutable)
            # translate node string names to indices
            self.index = {i:node for node, i in enumerate(self.nodes)}
            self.n = len(self.nodes)

        #self.functions = self._get_node_functions()
        #self.primes = self._get_primes()

        print('Imported Boolean graph with %i nodes'%self.n)


    def _squad_init(self, data):
        self.nodes = tuple(sorted(data.keys())) # tuple (i.e. immutable)
        self.n = len(self.nodes)
        self.index = {i:node for node, i in enumerate(self.nodes)}

        funcs = self._get_node_functions(data, self._squad_update_func)
        self.primes = QMC.functions2primes(funcs)


    def primes_to_matrices(self):
        '''
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Create logic matrices
        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to use e.g. dok_matrix
        '''
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        intgraph = InteractionGraphs.primes2igraph(self.primes)
        for node, d in intgraph.adjacency_iter(): #nx 1.11
            for other_node, sub_d in d.items():
                sign = next(iter(sub_d["sign"]))
                if sign == 1:
                    Act[self.index[other_node], self.index[node]] = 1
                elif sign == -1:
                    Inh[self.index[other_node], self.index[node]] = 1
                else:
                    print("Warning: Issue with edge: ", node, other_node)

        return ensure_ndarray(Act), ensure_ndarray(Inh)

    def interactions_to_matrices(self, data):
        '''
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Create logic matrices
        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to use e.g. dok_matrix
        '''
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        for node, d in data.items():
            for other_node, sign in d.items():
                if sign == "+":
                    Act[self.index[node], self.index[other_node]] = 1
                elif sign == "-":
                    Inh[self.index[node], self.index[other_node]] = 1
                else:
                    print("Warning: Issue with edge: ", node, other_node)
        return ensure_ndarray(Act), ensure_ndarray(Inh)



    def _squad_update_func(self, data, args, node):
        # TODO remove dep on act and inh

        Act, Inh = self.interactions_to_matrices(data)
        def func(*func_input):
            state = np.ones(self.n)
            for other_node, other_node_state in zip(args, func_input):
                state[self.index[other_node]] = other_node_state

            col_ones = np.ones((self.n))
            if Inh[self.index[node],:].dot(col_ones):
                inh = Inh[self.index[node],:].dot(state) > 0
            else:
                inh = 0
            if Act[self.index[node],:].dot(col_ones):
                act = Act[self.index[node],:].dot(state) > 0
            else:
                act = 1
            node_state = act * (1-inh)
            return node_state
        func.depends = args

        return func

    def _get_node_functions(self, data, update_func, default=1):
        funcs = {node:default for node in data.keys()}
        for node, d in data.items():
            args = list(d.keys())
            func = update_func(data, args, node)
            funcs[node] = func
        return funcs

    def generate_states(self, fixed={}):
        '''
        Generate all possible states of network of size n.
        fixed = a dictionary of nodes with fixed state
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
            for free_variable_state in product([0, 1], repeat=num_free_variables):
                this_state_array = state_array.copy()
                for node, state in zip(free_variables, free_variable_state):
                    this_state_array[self.index[node]] = state
                yield this_state_array

    def plot_transitions(self, fig, init=0):
        if isinstance(init, int) and len(init) == 1:
            init = [str(init)*self.n]

        stg = PyBoolNet.StateTransitionGraphs.primes2stg(self.primes, "synchronous", init)
        stg.graph["node"]["color"] = "cyan"
        stg.graph["node"]["height"] = 0.3
        stg.graph["node"]["width"] = 0.45


        PyBoolNet.StateTransitionGraphs.stg2image(stg, fig, LayoutEngine="dot")


    def bool_steady_states(self):
        steady_states = AspSolver.steady_states(self.primes)
        return steady_states

    def __print__():
        pass

    def __len__(self    ):
        return self.n

'''
    def _rename_variables(d, reverser, direction):

        if direction == 'forward':
            reformat_func = lambda x: "var-" + x.lower()
            new_dict = {node:reformat_func(node) for node in d.keys()}
            reverser = {reformat_func(node):node for node in d.keys()}
            return new_dict

        elif direction == 'backward':
            new_dict = {reverser[key]:value for key, value in d.items()}

        return new_dict

     def translate_primes(primes, reverser):
        new_primes = {}
        for node, d in primes.items():
            node = reverser[node]
            p_rev = []
            for sublist in d:
                sublist_rev = []
                for sub_d in sublist:
                    sub_d_rev = self._rename_variables(sub_d, reverser, 'backward')
                    sublist_rev.append(sub_d_rev)
                p_rev.append(sublist_rev)
            new_primes[node] = p_rev
        return new_primes
'''