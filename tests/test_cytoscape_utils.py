'''
Test booldog.utils.cytoscape_utils (with py4cytoscape's table reads mocked,
so no running Cytoscape instance is needed).
'''
import unittest
from unittest import mock

import pandas as pd

from booldog import BoolDogModel
from booldog.utils import cytoscape_utils

# X = (A & !B) | C, Y = X, as a logic circuit:
#   A -> and_0, B -> not_0 -> and_0, and_0 -> or_0, C -> or_0, or_0 -> X -> Y
BNET = "A, A\nB, B\nC, C\nX, (A & !B) | C\nY, X"


def _tables():
    g = BoolDogModel.from_bnet(BNET).to_networkx(as_logic_circuit=True)
    nodes = pd.DataFrame({"name": list(g.nodes),
                          "type": [g.nodes[n]["type"] for n in g.nodes]})
    edges = pd.DataFrame({"source": [u for u, _ in g.edges],
                          "target": [v for _, v in g.edges]})
    return nodes, edges


class TestExtractSubnetwork(unittest.TestCase):

    def setUp(self):
        self.nodes, self.edges = _tables()

        def get_table_columns(table="node", network=None):
            return self.edges if table == "edge" else self.nodes

        patcher = mock.patch.object(cytoscape_utils.p4c.tables,
                                    "get_table_columns",
                                    side_effect=get_table_columns)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_predecessors(self):
        self.assertEqual(
            cytoscape_utils.extract_subnetwork(1, ["X"]),
            {"X", "or_0", "and_0", "not_0", "A", "B", "C"})

    def test_successors_include_all_gate_inputs(self):
        '''Walking downstream from A reaches and_0 and or_0; their other
        inputs (not_0 -> B, and C) must be included too.'''
        self.assertEqual(
            cytoscape_utils.extract_subnetwork(1, ["A"],
                                               direction="successors"),
            {"A", "and_0", "not_0", "B", "or_0", "C", "X"})

    def test_verify_detects_incomplete_hyperedge(self):
        name_to_type = dict(zip(self.nodes["name"], self.nodes["type"]))
        missing = cytoscape_utils._find_incomplete_hyperedge(
            {"A", "and_0", "or_0", "C", "X"}, self.edges, name_to_type,
            ("species",))
        self.assertEqual(missing, ("and_0", "not_0"))

    def test_verify_accepts_complete_subgraph(self):
        name_to_type = dict(zip(self.nodes["name"], self.nodes["type"]))
        self.assertIsNone(
            cytoscape_utils._find_incomplete_hyperedge(
                {"X", "or_0", "and_0", "not_0", "A", "B", "C"}, self.edges,
                name_to_type, ("species",)))


if __name__ == '__main__':
    unittest.main()
