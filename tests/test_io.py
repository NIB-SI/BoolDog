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
from booldog.io import interaction_logic


class TestBoolDogModelFrom(unittest.TestCase):

    # ---- Boolean network formats ------------------------------------

    def test_from_bnet(self):
        bn = BoolDogModel.from_bnet(BooleanNetworkExamples.BNET)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    def test_from_bnet_header_without_space(self):
        # "targets,factors" (no space) is what pyboolnet's own bundled
        # binary and real-world bnet files (e.g. biodivine-boolean-models)
        # use; BoolNet's documented "targets, factors" (with space) must
        # also still work (see KNOWN_BUGS.md).
        bnet_no_space = BooleanNetworkExamples.BNET.replace(
            "targets, factors", "targets,factors")
        bn = BoolDogModel.from_bnet(bnet_no_space)
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

    def test_from_tabularqual_without_validation(self):
        bn = BoolDogModel.from_tabularqual(
            BooleanNetworkExamples.TABULARQUAL_FILE, validate=False)
        self.assertDictEqual(bn.primes, BooleanNetworkExamples.PRIMES)

    def test_from_tabularqual_passes_validate_on(self):
        from unittest import mock
        import booldog.io.tabularqual as tq

        for validate in (True, False):
            with mock.patch.object(tq, "convert_spreadsheet_to_sbml",
                                   wraps=tq.convert_spreadsheet_to_sbml) as convert:
                BoolDogModel.from_tabularqual(
                    BooleanNetworkExamples.TABULARQUAL_FILE, validate=validate)
            self.assertEqual(convert.call_args.kwargs["validate"], validate)

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

    def test_from_graphml_custom_logic(self):

        class ConstantLogic(interaction_logic.LogicBuilder):

            def build(self, node, regulators):
                return "1"

        bn = BoolDogModel.from_graphml(InteractionNetworkExamples.GRAPHML_FILE,
                                       edge_type_key="weight",
                                       logic=ConstantLogic())
        # every node in the example network is the target of some edge
        self.assertDictEqual(bn.primes,
                             {n: [[], [{}]] for n in bn.node_ids})

if __name__ == '__main__':
    unittest.main()
