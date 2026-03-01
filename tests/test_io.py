'''
Test Boolean network read/write functions.
'''
import sys
import unittest

import networkx as nx
import igraph as ig

from examples import BooleanNetworkExamples, InteractionNetworkExamples

# sys.path.append("../")
from booldog import BoolDogModel


class TestBoolDogModelFrom(unittest.TestCase):

    # ---- Boolean network formats ------------------------------------

    def test_from_bnet(self):
        bn = BoolDogModel.from_bnet(BooleanNetworkExamples.BNET)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    def test_from_primes(self):
        bn = BoolDogModel.from_primes(BooleanNetworkExamples.PRIMES)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    def test_from_sbmlqual(self):
        bn = BoolDogModel.from_sbmlqual(BooleanNetworkExamples.SBMLQUAL_FILE)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    def test_from_tabularqual(self):
        bn = BoolDogModel.from_tabularqual(
            BooleanNetworkExamples.TABULARQUAL_FILE)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    # ---- Interaction / graph formats --------------------------------

    def test_from_interactions(self):
        bn = BoolDogModel.from_interactions(
            InteractionNetworkExamples.INTERACTIONS,
            activator_symbol="+",
            inhibitor_symbol="-")

        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

    def test_from_sif(self):
        bn = BoolDogModel.from_sif(InteractionNetworkExamples.SIF_FILE,
                                   activator_symbol="1",
                                   inhibitor_symbol="-1")
        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

    def test_from_networkx(self):
        g = nx.DiGraph(InteractionNetworkExamples.DICT_OF_DICT)

        bn = BoolDogModel.from_networkx(g,
                                        activator_symbol="+",
                                        inhibitor_symbol="-",
                                        edge_type_key="interaction")
        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

    def test_from_igraph(self):

        g = ig.Graph.TupleList(InteractionNetworkExamples.INTERACTIONS,
                               directed=True,
                               edge_attrs=["interaction"])

        bn = BoolDogModel.from_igraph(g,
                                      activator_symbol="+",
                                      inhibitor_symbol="-",
                                      edge_type_key="interaction")
        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

    def test_from_graphml(self):
        bn = BoolDogModel.from_graphml(InteractionNetworkExamples.GRAPHML_FILE, edge_type_key="weight")
        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

    def test_from_graphml_yEd(self):
        bn = BoolDogModel.from_graphml(
            InteractionNetworkExamples.GRAPHML_YED_FILE,
            yEd_labels=True, yEd_arrow_head=True)
        self.assertDictEqual(bn.primes,
                             InteractionNetworkExamples.PRIMES_SQUAD)

if __name__ == '__main__':
    unittest.main()
