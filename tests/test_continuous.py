'''
Test continuous/semi-quantitative simulation (booldog.continuous).
'''
import itertools
import os
import tempfile
import unittest

import numpy as np

from examples import BooleanNetworkExamples
from booldog import BoolDogModel
from booldog.continuous.ode_factory import ode_factory, BooleCubeODE, SquadODE


def _model():
    return BoolDogModel.from_bnet(BooleanNetworkExamples.BNET)


def _regulator_edge_case_model():
    return BoolDogModel.from_bnet(BooleanNetworkExamples.BNET_REGULATOR_EDGE_CASES)


def _boolean_function_value(primes_for_node, state):
    '''Evaluate a node's Boolean function directly from its prime
    implicants (DNF over the "on" primes), independent of any ODE
    machinery - the ground truth B1 is interpolating.'''
    return float(any(
        all(state[var] == val for var, val in prime.items())
        for prime in primes_for_node[1]))


class TestOdeFactory(unittest.TestCase):

    def test_unknown_transform_raises(self):
        '''An invalid transform name raises ValueError.'''
        with self.assertRaises(ValueError):
            ode_factory(_model(), "not-a-real-transform")

    def test_transform_is_case_insensitive(self):
        '''Transform names are matched case-insensitively.'''
        ode = ode_factory(_model(), "BoolECube")
        self.assertIsInstance(ode, BooleCubeODE)

    def test_boolecube_family_dispatch(self):
        '''boolecube/hillcube/normalisedhillcube all route to BooleCubeODE.'''
        for transform in ["boolecube", "hillcube", "normalisedhillcube"]:
            ode = ode_factory(_model(), transform)
            self.assertIsInstance(ode, BooleCubeODE)
            self.assertEqual(ode.transform, transform)

    def test_squad_dispatch(self):
        '''squad routes to SquadODE.'''
        ode = ode_factory(_model(), "squad")
        self.assertIsInstance(ode, SquadODE)


class TestODEBase(unittest.TestCase):

    def test_event_function(self):
        '''event_function returns t - event_t, the trigger solve_ivp uses to
        stop exactly at a perturbation time.'''
        ode = ode_factory(_model(), "boolecube")
        self.assertEqual(ode.event_function(5, None, 2), 3)


class TestBooleCubeODE(unittest.TestCase):

    def setUp(self):
        self.model = _model()

    def test_params_as_scalars(self):
        '''Scalar tau/n/k broadcast to per-node arrays, and boolecube wires
        up identity as the transform function.'''
        ode = ode_factory(self.model, "boolecube", tau=1, n=3, k=0.5)
        self.assertEqual(ode.n, len(self.model))
        np.testing.assert_array_equal(ode.param_tau, np.ones(ode.n))
        np.testing.assert_array_equal(ode.param_n, np.full(ode.n, 3))
        np.testing.assert_array_equal(ode.param_k, np.full(ode.n, 0.5))
        self.assertEqual(ode.transform_function, ode.identity)

    def test_hillcube_uses_hill_transform(self):
        '''hillcube wires up hill as the transform function.'''
        ode = ode_factory(self.model, "hillcube")
        self.assertEqual(ode.transform_function, ode.hill)

    def test_normalisedhillcube_uses_normalised_hill_transform(self):
        '''normalisedhillcube wires up normalised_hill as the transform
        function.'''
        ode = ode_factory(self.model, "normalisedhillcube")
        self.assertEqual(ode.transform_function, ode.normalised_hill)

    def test_identity(self):
        '''identity is a no-op.'''
        ode = ode_factory(self.model, "boolecube")
        x = np.array([0.1, 0.9] + [0.5] * (ode.n - 2))
        np.testing.assert_array_equal(ode.identity(x), x)

    def test_hill_function(self):
        '''hill matches the hand-computed Hill function.'''
        ode = ode_factory(self.model, "boolecube", n=2, k=0.5)
        x = np.full(ode.n, 0.5)
        expected = x**2 / (x**2 + 0.5**2)
        np.testing.assert_allclose(ode.hill(x), expected)

    def test_normalised_hill_function(self):
        '''normalised_hill is hill(x) normalised by hill(1).'''
        ode = ode_factory(self.model, "boolecube", n=2, k=0.5)
        x = np.full(ode.n, 0.5)
        np.testing.assert_allclose(ode.normalised_hill(x),
                                   ode.hill(x) / ode.hill(np.ones(ode.n)))

    def test_B1_matches_boolean_function_at_every_corner(self):
        '''B1 is a multilinear interpolation of the network's prime
        implicants, so at every corner of the Boolean state space (i.e.
        x in {0,1}^n) it must exactly reproduce the underlying Boolean
        function - checked here against primes directly, independent of
        any ODE/integration machinery, so this is scipy-free: it tests
        that the ODE was *built* correctly, not how it integrates.'''
        ode = ode_factory(self.model, "boolecube")
        for corner in itertools.product([0, 1], repeat=ode.n):
            state = {node: corner[self.model.index[node]]
                     for node in self.model.node_ids}
            b1 = ode.B1(np.array(corner, dtype=float))
            for node in self.model.node_ids:
                expected = _boolean_function_value(self.model.primes[node], state)
                with self.subTest(corner=corner, node=node):
                    self.assertEqual(b1[self.model.index[node]], expected)

    def test_dxdt_shape(self):
        '''dxdt returns one derivative value per node.'''
        ode = ode_factory(self.model, "boolecube")
        x = np.full(ode.n, 0.5)
        d = ode.dxdt(0, x.copy())
        self.assertEqual(d.shape, (ode.n,))

    def test_dxdt_clips_out_of_unit_range(self):
        '''States outside [0, 1] are clipped internally rather than
        producing inf/nan.'''
        ode = ode_factory(self.model, "boolecube")
        x = np.array([-0.5, 1.5] + [0.5] * (ode.n - 2))
        d = ode.dxdt(0, x.copy())
        self.assertEqual(d.shape, (ode.n,))
        self.assertTrue(np.all(np.isfinite(d)))

    def test_zero_tau_freezes_node(self):
        '''A node with tau=0 always has derivative 0.'''
        tau = {name: 1 for name in self.model.node_ids}
        first_node = self.model.node_ids[0]
        tau[first_node] = 0
        ode = ode_factory(self.model, "boolecube", tau=tau)
        d = ode.dxdt(0, np.full(ode.n, 0.7))
        self.assertEqual(d[self.model.index[first_node]], 0)

    def test_update_sets_off_nodes(self):
        '''update(off_nodes=...) freezes a node's derivative to 0, the
        mechanism continuous_simulation uses for timed perturbations.'''
        ode = ode_factory(self.model, "boolecube")
        idx = self.model.index[self.model.node_ids[0]]
        ode.update(off_nodes={idx})
        d = ode.dxdt(0, np.full(ode.n, 0.7))
        self.assertEqual(d[idx], 0)


class TestSquadODE(unittest.TestCase):

    def setUp(self):
        self.model = _model()
        self.ode = ode_factory(self.model, "squad", gamma=1, h=10)

    def test_params_as_scalars(self):
        '''Scalar gamma/h broadcast to per-node arrays.'''
        np.testing.assert_array_equal(self.ode.param_gamma, np.ones(self.ode.n))
        np.testing.assert_array_equal(self.ode.param_h, np.full(self.ode.n, 10))

    def test_activation_inhibition_matrix_values(self):
        '''activations/inhibitions correctly classify each node's regulators:
        node_C is activated by node_A and node_B only (its rule is
        "A & B"); node_E is activated by node_D and inhibited by node_C
        (its rule is "D | !C").'''
        act, inh = self.ode.activations, self.ode.inhibitions
        idx = self.model.index

        expected_act_c = np.zeros(self.ode.n)
        expected_act_c[idx["node_A"]] = 1
        expected_act_c[idx["node_B"]] = 1
        np.testing.assert_array_equal(act[idx["node_C"]], expected_act_c)
        np.testing.assert_array_equal(inh[idx["node_C"]], np.zeros(self.ode.n))

        expected_act_e = np.zeros(self.ode.n)
        expected_act_e[idx["node_D"]] = 1
        expected_inh_e = np.zeros(self.ode.n)
        expected_inh_e[idx["node_C"]] = 1
        np.testing.assert_array_equal(act[idx["node_E"]], expected_act_e)
        np.testing.assert_array_equal(inh[idx["node_E"]], expected_inh_e)

    def test_prefactors_match_activator_inhibitor_counts(self):
        '''_a1/_b1 are (1 + count) / count for a node's number of
        activators/inhibitors (see the comment in SquadODE.__init__ on why
        this differs from the literal, seemingly erratum'd, 2007 paper).'''
        idx = self.model.index
        # node_C: 2 activators (A, B), 0 inhibitors
        self.assertAlmostEqual(self.ode._a1[idx["node_C"]], 1.5)
        # node_E: 1 activator (D), 1 inhibitor (C)
        self.assertAlmostEqual(self.ode._a1[idx["node_E"]], 2.0)
        self.assertAlmostEqual(self.ode._b1[idx["node_E"]], 2.0)

    def test_omega_bounded_in_unit_interval(self):
        '''_omega (Di Cara et al. 2007, eq. 2) stays within [0, 1].'''
        w = self.ode._omega(np.full(self.ode.n, 0.5))
        self.assertTrue(np.all(w >= 0) and np.all(w <= 1))

    def test_omega_activators_only(self):
        '''node_C (2 activators, 0 inhibitors) matches the "only
        activators" case of eq. 2 exactly: a1 * (a_x / (1 + a_x)).'''
        idx = self.model.index
        x = np.zeros(self.ode.n)
        x[idx["node_A"]] = 1
        x[idx["node_B"]] = 1
        w = self.ode._omega(x)
        self.assertAlmostEqual(w[idx["node_C"]], 1.0)

        x[idx["node_B"]] = 0
        w = self.ode._omega(x)
        # only one of two activators on: a1 * (1/2) = 1.5 * 0.5 = 0.75.
        # Not 0, unlike the true AND - this is the known SQUAD limitation
        # discussed earlier (it aggregates activators additively, it
        # cannot express "requires all activators").
        self.assertAlmostEqual(w[idx["node_C"]], 0.75)

    def test_omega_inhibitors_only(self):
        '''C (0 activators, 1 inhibitor: A) matches the "only inhibitors"
        case of eq. 2 exactly: 1 - b1 * (b_x / (1 + b_x)).'''
        model = _regulator_edge_case_model()
        ode = ode_factory(model, "squad", gamma=1, h=10)
        idx = model.index
        x = np.zeros(ode.n)
        x[idx["A"]] = 0.5
        w = ode._omega(x)
        # b1 = (1+1)/1 = 2, b_x = 0.5: 1 - 2*(0.5/1.5) = 1/3.
        self.assertAlmostEqual(w[idx["C"]], 1 / 3)

    def test_omega_activators_and_inhibitors(self):
        '''node_E (1 activator: D, 1 inhibitor: C) matches the "both" case
        of eq. 2 exactly: a1*(a_x/(1+a_x)) * (1 - b1*(b_x/(1+b_x))).'''
        idx = self.model.index

        def omega_e(d_value, c_value):
            x = np.zeros(self.ode.n)
            x[idx["node_D"]] = d_value
            x[idx["node_C"]] = c_value
            return self.ode._omega(x)[idx["node_E"]]

        # fully activated, no inhibition: a1*(1/2)=1 times (1-0)=1 -> 1.0
        self.assertAlmostEqual(omega_e(1, 0), 1.0)
        # fully activated AND fully inhibited: activation term is 1, but
        # b1*(1/2)=1 makes the inhibition term (1-1)=0 -> cancels to 0.
        self.assertAlmostEqual(omega_e(1, 1), 0.0)
        # activator off: activation term is 0 regardless of inhibitor.
        self.assertAlmostEqual(omega_e(0, 0), 0.0)
        self.assertAlmostEqual(omega_e(0, 1), 0.0)

    def test_omega_fully_unregulated_node_is_zero(self):
        '''D has no regulators at all (0 activators, 0 inhibitors) - the
        explicit override in _omega forces this to exactly 0, rather than
        whatever the masked activator/inhibitor terms would otherwise give.'''
        model = _regulator_edge_case_model()
        ode = ode_factory(model, "squad", gamma=1, h=10)
        idx = model.index
        for value in [0, 0.3, 0.5, 1]:
            with self.subTest(value=value):
                x = np.full(ode.n, value)
                w = ode._omega(x)
                self.assertEqual(w[idx["D"]], 0)

    def test_dxdt_transform_boundary_normalisation(self):
        '''_dxdt_transform's sigmoid term is exactly 0 at w=0 and exactly 1
        at w=1 (isolated here via x=0, so the -gamma*x decay term drops
        out) - the boundary-exact property that distinguishes this from
        the literal (non-boundary-exact) 2007 formula, see the comment on
        _dxdt_transform.'''
        zeros = np.zeros(self.ode.n)
        np.testing.assert_allclose(
            self.ode._dxdt_transform(zeros, np.zeros(self.ode.n)), 0, atol=1e-10)
        np.testing.assert_allclose(
            self.ode._dxdt_transform(zeros, np.ones(self.ode.n)), 1, atol=1e-10)

    def test_dxdt_shape(self):
        '''dxdt returns one derivative value per node.'''
        d = self.ode.dxdt(0, np.full(self.ode.n, 0.5))
        self.assertEqual(d.shape, (self.ode.n,))

    def test_update_sets_off_nodes(self):
        '''update(off_nodes=...) freezes a node's derivative to 0.'''
        idx = 0
        self.ode.update(off_nodes=[idx])
        d = self.ode.dxdt(0, np.full(self.ode.n, 0.5))
        self.assertEqual(d[idx], 0)


class TestTransformBoolToContinuous(unittest.TestCase):

    def test_default_transform(self):
        '''Default transform is normalisedhillcube.'''
        ode = _model().transform_bool_to_continuous()
        self.assertIsInstance(ode, BooleCubeODE)
        self.assertEqual(ode.transform, "normalisedhillcube")

    def test_explicit_transform_and_kwargs(self):
        '''An explicit transform and its kwargs (e.g. gamma for squad) are
        passed through to the ODE class.'''
        ode = _model().transform_bool_to_continuous(transform="squad", gamma=2)
        self.assertIsInstance(ode, SquadODE)
        np.testing.assert_array_equal(ode.param_gamma, np.full(ode.n, 2))


class TestContinuousSimulation(unittest.TestCase):

    def setUp(self):
        self.model = _model()

    def test_basic_run_shapes(self):
        '''A basic run's result spans t_min to t_max with one column of y
        per node.'''
        result = self.model.continuous_simulation(
            t_min=0, t_max=2, initial_state=0.5, transform="boolecube")
        self.assertEqual(result.y.shape[1], len(self.model))
        self.assertEqual(result.y.shape[0], result.t.shape[0])
        self.assertAlmostEqual(result.t[0], 0)
        self.assertAlmostEqual(result.t[-1], 2)

    def test_initial_state_dict_with_default(self):
        '''A dict initial_state sets named nodes explicitly and the rest via
        its "default" key.'''
        first_node = self.model.node_ids[0]
        result = self.model.continuous_simulation(
            t_min=0, t_max=1,
            initial_state={first_node: 0.9, "default": 0.1},
            transform="boolecube")
        self.assertAlmostEqual(result.y[0, self.model.index[first_node]], 0.9)

    def test_node_events_point_perturbation_sets_exact_value(self):
        '''A one-off node_event (no duration) sets the node to exactly
        "value" at that instant, then lets it evolve freely again
        (simulation still reaches t_max).'''
        node = self.model.node_ids[0]
        idx = self.model.index[node]
        result = self.model.continuous_simulation(
            t_min=0, t_max=3, initial_state=0,
            node_events=[{"time": 1, "node": node, "value": 0.6}],
            transform="boolecube")
        self.assertGreaterEqual(result.t[-1], 3 - 1e-6)
        # solve_ivp restarts exactly at t=1, so the boundary point is
        # duplicated: the last of the ties is the post-perturbation value.
        at_event = np.where(np.isclose(result.t, 1))[0][-1]
        self.assertAlmostEqual(result.y[at_event, idx], 0.6)

    def test_node_events_with_duration_holds_exact_value_then_releases(self):
        '''A node_event with a duration freezes the node at exactly "value"
        for the whole window (derivative forced to 0), then releases it to
        evolve freely again once the window ends.'''
        node = self.model.node_ids[0]
        idx = self.model.index[node]
        result = self.model.continuous_simulation(
            t_min=0, t_max=3, initial_state=0,
            node_events=[{"time": 1, "node": node, "value": 0.75, "duration": 1}],
            transform="boolecube")
        self.assertGreaterEqual(result.t[-1], 3 - 1e-6)
        # t=1 is a duplicated boundary point (pre- and post-event value);
        # strictly-inside points are unambiguous.
        in_window = (result.t > 1) & (result.t <= 2)
        self.assertTrue(in_window.any())
        np.testing.assert_allclose(result.y[in_window, idx], 0.75, atol=1e-9)
        # released afterward: no longer pinned exactly to the forced value.
        self.assertFalse(np.isclose(result.y[-1, idx], 0.75, atol=1e-9))

    def test_single_node_event_as_dict_not_list(self):
        '''A single node_event may be passed as a bare dict instead of a
        list of dicts.'''
        node = self.model.node_ids[0]
        result = self.model.continuous_simulation(
            t_min=0, t_max=1, initial_state=0,
            node_events={"time": 0.5, "node": node, "value": 1},
            transform="boolecube")
        self.assertGreaterEqual(result.t[-1], 1 - 1e-6)

    def test_converges_to_expected_boolean_steady_state(self):
        '''The continuous relaxation actually approximates the underlying
        Boolean function: node_C's rule is "node_A & node_B" (see
        BooleanNetworkExamples.PRIMES), so holding A and B fixed (tau=0, a
        steep hillcube) should drive C to ~1 only when both are 1, and to
        ~0 otherwise - not just "doesn't crash".'''
        tau = {name: 1 for name in self.model.node_ids}
        tau["node_A"] = 0
        tau["node_B"] = 0
        idx_c = self.model.index["node_C"]

        for a_value, b_value, expected in [(1, 1, 1), (0, 1, 0),
                                           (1, 0, 0), (0, 0, 0)]:
            with self.subTest(node_A=a_value, node_B=b_value):
                result = self.model.continuous_simulation(
                    t_min=0, t_max=5,
                    initial_state={"node_A": a_value, "node_B": b_value, "default": 0},
                    transform="hillcube", tau=tau, n=10, k=0.5)
                self.assertAlmostEqual(result.y[-1, idx_c], expected, places=1)

    def test_reuses_provided_ode_system(self):
        '''Passing ode_system= skips building a new ODE system and uses the
        one given.'''
        ode = self.model.transform_bool_to_continuous(transform="boolecube")
        result = self.model.continuous_simulation(
            t_min=0, t_max=1, initial_state=0, ode_system=ode)
        self.assertIs(result.ode_system, ode)

    def test_export(self):
        '''export() writes the transform name and node list to the output
        file.'''
        result = self.model.continuous_simulation(
            t_min=0, t_max=1, initial_state=0.2, transform="boolecube")
        with tempfile.TemporaryDirectory() as tmp_dir:
            outfile = os.path.join(tmp_dir, "out.tsv")
            result.export(outfile)
            with open(outfile, encoding="utf-8") as f:
                content = f.read()
        self.assertIn("#transform\tboolecube", content)
        self.assertIn("#nodelist\t", content)


if __name__ == '__main__':
    unittest.main()
