import os
import sys
import unittest

import networkx as nx
import igraph as ig

from booldog.io import read

example_primes_true = {
    'node_A': [
        [{'node_B': 0, 'node_E': 0}, {'node_B': 1, 'node_E': 1}],
        [{'node_B': 0, 'node_E': 1}, {'node_B': 1, 'node_E': 0}]
    ],
    'node_B': [
        [{'node_B': 0}],
        [{'node_B': 1}]
    ],
    'node_C': [
        [{'node_B': 0}, {'node_A': 0}],
        [{'node_A': 1, 'node_B': 1}]], 'node_D': [[{'node_D': 0}], [{'node_D': 1}]
    ],
    'node_E': [
        [{'node_C': 1, 'node_D': 0}],
        [{'node_D': 1}, {'node_C': 0}]
    ]
}
example_primes_squad = {
    'node_A': [
        [{}],
        []
    ],
    'node_B': [
        [{'node_B': 0}],
        [{'node_B': 1}]
    ],
    'node_C': [
        [{'node_A': 0, 'node_B': 0}],
        [{'node_B': 1}, {'node_A': 1}]
    ],
    'node_D': [
        [{'node_D': 0}], [{'node_D': 1}]
    ],
    'node_E': [
        [{'node_D': 0}, {'node_C': 1}],
        [{'node_C': 0, 'node_D': 1}]
    ]
}
example_interactions = {
    'node_A': {
        'node_C': '+'
        },
    'node_B': {
        'node_A': '-',
        'node_B': '+',
        'node_C': '+'
    },
    'node_C': {
        'node_E': '-'
    },
    'node_D': {
        'node_D': '+',
        'node_E': '+'
    },
    'node_E': {
        'node_A': '-'
    }
}
example_edgelist = [(s, t, example_interactions[s][t]) for s in example_interactions for t in example_interactions[s]]
example_dict_of_dict = {s:{t:{"interaction": example_interactions[s][t]}} for s in example_interactions for t in example_interactions[s]}

class Test(unittest.TestCase):

    def test_read_bnet(self):
        G = read("data/example.bnet", "bnet")
        self.assertDictEqual(G.primes, example_primes_true)

    def test_read_primes(self):
        G = read("data/example.primes.json", "primes")
        self.assertDictEqual(G.primes, example_primes_true)

    def test_read_sbmlqual(self):
        G = read("data/example.xml", "sbml-qual")
        self.assertDictEqual(G.primes, example_primes_true)

    ##################################################
    # Uses SQUAD to turn into interactions Boolean Graph
    ##################################################

    def test_read_sif(self):
        G = read("data/example.sif", "sif", activator_value="1", inhibitor_value="-1")
        self.assertDictEqual(G.primes, example_primes_squad)

    def test_read_interactions(self):
        G = read(example_interactions, "interactions", activator_value="+", inhibitor_value="-")
        self.assertDictEqual(G.primes, example_primes_squad)

    def test_read_graphml(self):
        G = read("data/example.graphml", "graphml")
        self.assertDictEqual(G.primes, example_primes_squad)

    def test_read_graphml_yEd(self):
        G = read("data/example.yEd.graphml", "graphml", yEd=True, use_labels=False)
        self.assertDictEqual(G.primes, example_primes_squad)

    def test_read_igraph(self):
        graph = ig.Graph.TupleList(example_edgelist, directed=True, edge_attrs=["weight"])
        G = read(graph, "igraph", activator_value="+", inhibitor_value="-")
        self.assertDictEqual(G.primes, example_primes_squad)

    def test_read_networkx(self):
        graph = nx.DiGraph(example_dict_of_dict)
        G = read(graph, "networkx", activator_value="+", inhibitor_value="-", edge_type_key="interaction")
        self.assertDictEqual(G.primes, example_primes_squad)

if __name__ == '__main__':
    unittest.main()
