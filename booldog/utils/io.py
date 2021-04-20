import os
import sys

import xmltodict
import igraph as ig

import numpy as np
from .utils import *


from PyBoolNet import FileExchange
from pathlib import Path

#from PyBoolNet import QuineMcCluskey as QMC
from . import QMC

format_classes = {'primes':0,
                  'bnet':0,
                  'interactions':0,
                  'graphml':0
                  }

def read_boolean_graph(graph, in_format, default=1):


    if not (in_format in format_classes):
       raise ValueError(f"'in_format' arguments must be in "\
                        f"{list(format_classes.keys())}")

    if not isinstance(graph, (dict, str, Path)):
       raise TypeError("'graph' argument must be a "\
                       "dictionary, string or pathlib.Path object")

    # if isinstance(graph, (str, Path)):
    #     file_path = Path(graph)
    #     if not file_path.exists():
    #        raise FileNotFoundError(
    #             f"'graph' argument {graph} is a path or string, "\
    #             f"but does not exist.")

    if (in_format == 'primes') and isinstance(graph, dict):
        # nothing to do
        primes = graph

    elif (in_format == 'interactions') and isinstance(graph, dict):
        # TODO: This is ugly, try some inheritance?
        g = SquadInteractions(graph, default=default)
        primes, nodes, index, n = g.primes, g.nodes, g.index, g.n

    elif (in_format == 'primes') and (os.path.isfile(graph)):
        # primes from a file
        primes = FileExchange.read_primes(graph)

    elif (in_format == 'bnet') and \
         (isinstance(graph, str) or os.path.isfile(graph)):
        # primes from a file
        primes = FileExchange.bnet2primes(graph)

    elif (in_format == 'graphml') and (os.path.isfile(graph)):
        primes = import_graphml(path, **kwargs)

    else:
        raise NotImplementedError(f"import of {in_format} from {type(graph)} "\
                                  f"is not implemented. ")

    if not (in_format == 'interactions'):
        nodes = tuple(sorted(primes.keys())) # tuple (i.e. immutable)
        index = {i:node for node, i in enumerate(nodes)}
        n = len(nodes)

        #self.functions = self._get_node_functions()
        #self.primes = self._get_primes()

    return primes, nodes, index, n


def write_boolean_graph(graph, out_format, outfile):

    if not  (out_format in formats):
       raise ValueError(f"'out_format' argument must be in "\
                        f"{list(format_classes.keys())}")

    outfile = Path(outfile)
    file_writable(outfile)

    #TODO


##############################
#     SQUAD/INTERACTIONS     #
##############################

class SquadInteractions:

    def __init__(self, data, default=1):
        self.nodes = tuple(sorted(data.keys())) # tuple (i.e. immutable)
        self.n = len(self.nodes)
        self.index = {i:node for node, i in enumerate(self.nodes)}

        funcs = self.squad_update_funcs(data, default=default)
        self.primes = QMC.functions2primes(funcs)


    def interactions_to_matrices(self, data):
        '''
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Create logic matrices
        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to
        use e.g. dok_matrix
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

    def squad_update_funcs(self, data, default=1):
        # TODO remove dep on act and inh
        funcs = {node:default for node in data.keys()}

        Act, Inh = self.interactions_to_matrices(data)

        for node, d in data.items():
            args = list(d.keys())

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

            funcs[node] = func

        return funcs




##############################
#           OTHER            #
##############################

def import_graphml(path, inhibitor_symbol="white_diamond", activator_symbol="standard"):
    # load graph
    g = ig.Graph.Read_GraphML(path)

    # add edge attributes (i.e. activator or inhibitor)
    with open(path) as f:
        d = xmltodict.parse(f.read())

    D = {}
    N = {}

    for v in d["graphml"]["graph"]["node"]:
        N[v["@id"]] = v["data"]["y:ShapeNode"]["y:NodeLabel"]["#text"]

    for e in d["graphml"]["graph"]["edge"]:
        symbol = e["data"]["y:PolyLineEdge"]['y:Arrows']['@target']
        if symbol == activator_symbol:
            D[e["@id"]] = "+"
        elif symbol == inhibitor_symbol:
            D[e["@id"]] = "-"
        else:
            print("Issue with edge ", e["@id"], "symbol not activator or inhibitor", symbol)

    for e in g.es():
        e["type"] = D[e["id"]]

    for v in g.vs():
        v["id"] = N[v["id"]]

    g_dict = {node["id"]:{} for node in g.vs()}
    for e in g.es():
        source = g.vs()[e.source]["id"]
        target = g.vs()[e.target]["id"]
        sign = e["type"]
        g_dict[target][source] = sign

    return g_dict
