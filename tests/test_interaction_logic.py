''' Tests for booldog.io.interaction_logic
'''
import unittest

from booldog.io import interaction_logic

LOGGER_NAME = "booldog.io.interaction_logic"


class TestSquadLogic(unittest.TestCase):

    def setUp(self):
        self.logic = interaction_logic.SquadLogic()

    def test_no_regulators(self):
        self.assertEqual(self.logic.build("X", {}), "0")

    def test_activators_only(self):
        self.assertEqual(self.logic.build("X", {"A": 1, "B": 1}), "A | B")

    def test_inhibitors_only(self):
        self.assertEqual(self.logic.build("X", {"C": -1, "D": -1}),
                         "!C & !D")

    def test_activators_and_inhibitors(self):
        self.assertEqual(
            self.logic.build("X", {"A": 1, "C": -1, "B": 1, "D": -1}),
            "(A | B) & (!C & !D)")


class TestInteractions2Rules(unittest.TestCase):

    def test_unrecognised_sign_dropped_with_warning(self):
        with self.assertLogs(LOGGER_NAME, level="WARNING") as cm:
            rules = interaction_logic.interactions2rules([("A", "X", 1),
                                                          ("B", "X", 0)])
        self.assertEqual(rules, {"X": "A"})
        self.assertEqual(len(cm.records), 1)
        self.assertIn("unrecognized sign 0", cm.output[0])

    def test_duplicate_conflicting_last_kept_with_warning(self):
        with self.assertLogs(LOGGER_NAME, level="WARNING") as cm:
            rules = interaction_logic.interactions2rules([("A", "X", 1),
                                                          ("A", "X", -1)])
        self.assertEqual(rules, {"X": "!A"})
        self.assertEqual(len(cm.records), 1)
        self.assertIn("Conflicting signs", cm.output[0])

    def test_duplicate_same_sign_last_kept_with_warning(self):
        with self.assertLogs(LOGGER_NAME, level="WARNING") as cm:
            rules = interaction_logic.interactions2rules([("A", "X", -1),
                                                          ("A", "X", -1)])
        self.assertEqual(rules, {"X": "!A"})
        self.assertEqual(len(cm.records), 1)
        self.assertIn("Repeated edge", cm.output[0])
        self.assertNotIn("Conflicting", cm.output[0])

    def test_custom_logic(self):

        class AndLogic(interaction_logic.LogicBuilder):

            def build(self, node, regulators):
                return " & ".join(regulators)

        rules = interaction_logic.interactions2rules(
            [("A", "X", 1), ("B", "X", 1)], logic=AndLogic())
        self.assertEqual(rules, {"X": "A & B"})


if __name__ == '__main__':
    unittest.main()
