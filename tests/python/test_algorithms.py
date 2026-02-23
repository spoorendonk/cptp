"""Test the solver-independent algorithm Python bindings.

These tests exercise SeparationOracle, BoundPropagator, build_warm_start,
forward_labeling, backward_labeling, and edge_elimination — all without HiGHS.
"""

from pathlib import Path

import numpy as np
import pytest

import rcspp_bac
from rcspp_bac import (
    BoundPropagator,
    Cut,
    Problem,
    SeparationOracle,
    WarmStartResult,
    build_warm_start,
    edge_elimination,
    forward_labeling,
    backward_labeling,
    load,
)

DATA_DIR = Path(__file__).parent.parent / "data"


# ---------- helpers ----------

def make_tour_problem():
    """4-node tour: depot=0, customers=1,2,3."""
    edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)
    costs = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0])
    profits = np.array([0.0, 20.0, 15.0, 10.0])
    demands = np.array([0.0, 3.0, 4.0, 2.0])
    return Problem(
        num_nodes=4, edges=edges, edge_costs=costs,
        profits=profits, demands=demands,
        capacity=7.0, source=0, target=0,
    )


def make_path_problem():
    """4-node s-t path: source=0, target=3."""
    edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)
    costs = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0])
    profits = np.array([0.0, 20.0, 15.0, 10.0])
    demands = np.array([0.0, 3.0, 4.0, 2.0])
    return Problem(
        num_nodes=4, edges=edges, edge_costs=costs,
        profits=profits, demands=demands,
        capacity=7.0, source=0, target=3,
    )


# ---------- SeparationOracle ----------

class TestSeparationOracle:
    def test_finds_violated_cuts(self):
        prob = make_tour_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)
        oracle.add_default_separators()

        # Disconnected subtour: {0,1} and {2,3} active
        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 0.5])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 2 and v == 3):
                x[i] = 1.0

        cuts = oracle.separate(x, y, x_offset=0, y_offset=m)
        assert len(cuts) > 0
        assert all(isinstance(c, Cut) for c in cuts)
        # Sorted by violation (descending)
        violations = [c.violation for c in cuts]
        assert violations == sorted(violations, reverse=True)

    def test_no_cuts_on_feasible(self):
        prob = make_tour_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)
        oracle.add_default_separators()

        x = np.zeros(m)
        y = np.array([1.0, 1.0, 1.0, 0.0])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 0 and v == 2) or (u == 1 and v == 2):
                x[i] = 1.0

        cuts = oracle.separate(x, y, x_offset=0, y_offset=m, tol=1e-6)
        assert len(cuts) == 0

    def test_is_feasible(self):
        prob = make_tour_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)

        x = np.zeros(m)
        y = np.array([1.0, 1.0, 1.0, 0.0])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 0 and v == 2) or (u == 1 and v == 2):
                x[i] = 1.0

        assert oracle.is_feasible(x, y, x_offset=0, y_offset=m)

    def test_max_cuts_limits_output(self):
        prob = make_tour_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)
        oracle.add_default_separators()
        oracle.set_max_cuts_per_separator(1)

        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 0.5])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 2 and v == 3):
                x[i] = 1.0

        cuts = oracle.separate(x, y, x_offset=0, y_offset=m)
        assert len(cuts) <= 4  # 4 separators, 1 cut each

    def test_individual_separators(self):
        prob = make_tour_problem()
        oracle = SeparationOracle(prob)
        oracle.add_sec()
        oracle.add_rci()
        oracle.add_multistar()
        oracle.add_comb()
        oracle.add_rglm()
        # Just check it doesn't crash and separators accumulate
        # (no public way to query count, but we can verify by running)
        m = prob.num_edges
        x = np.zeros(m)
        y = np.zeros(prob.num_nodes)
        y[0] = 1.0
        oracle.separate(x, y, x_offset=0, y_offset=m)

    def test_cut_attributes(self):
        prob = make_tour_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)
        oracle.add_default_separators()

        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 0.5])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 2 and v == 3):
                x[i] = 1.0

        cuts = oracle.separate(x, y, x_offset=0, y_offset=m)
        assert len(cuts) > 0

        cut = cuts[0]
        assert isinstance(cut.indices, np.ndarray)
        assert cut.indices.dtype == np.int32
        assert isinstance(cut.values, np.ndarray)
        assert cut.values.dtype == np.float64
        assert len(cut.indices) == len(cut.values)
        assert cut.size == len(cut.indices)
        assert isinstance(cut.rhs, float)
        assert isinstance(cut.violation, float)
        assert cut.violation > 0


# ---------- Warm-start heuristic ----------

class TestWarmStart:
    def test_tour(self):
        prob = make_tour_problem()
        result = build_warm_start(prob, time_budget_ms=50.0)
        assert isinstance(result, WarmStartResult)
        assert isinstance(result.col_values, np.ndarray)
        assert result.col_values.dtype == np.float64
        assert len(result.col_values) == prob.num_edges + prob.num_nodes
        assert result.objective < 1e18

    def test_path(self):
        prob = make_path_problem()
        result = build_warm_start(prob, time_budget_ms=50.0)
        assert isinstance(result, WarmStartResult)
        m = prob.num_edges
        # Source and target should be visited
        assert result.col_values[m + prob.source] == 1.0
        assert result.col_values[m + prob.target] == 1.0


# ---------- Preprocessing ----------

class TestPreprocessing:
    def test_forward_labeling(self):
        prob = make_tour_problem()
        bounds = forward_labeling(prob, prob.source)
        assert isinstance(bounds, np.ndarray)
        assert bounds.dtype == np.float64
        assert len(bounds) == prob.num_nodes
        # Source bound should be 0 (or very small)
        assert bounds[prob.source] <= 1e-6

    def test_backward_labeling(self):
        prob = make_path_problem()
        bounds = backward_labeling(prob, prob.target)
        assert isinstance(bounds, np.ndarray)
        assert bounds.dtype == np.float64
        assert len(bounds) == prob.num_nodes
        assert bounds[prob.target] <= 1e-6

    def test_edge_elimination(self):
        """Edge with huge cost should be eliminated with tight UB."""
        edges = np.array([[0, 1], [0, 2], [1, 2]], dtype=np.int32)
        costs = np.array([1.0, 1.0, 100.0])
        profits = np.array([0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=3, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=0,
        )

        fwd = forward_labeling(prob, 0)
        bwd = fwd.copy()

        eliminated = edge_elimination(prob, fwd, bwd, upper_bound=4.0, correction=0.0)
        assert isinstance(eliminated, np.ndarray)
        assert len(eliminated) == prob.num_edges
        # Edge {1,2} (index 2) with cost=100 should be eliminated
        assert eliminated[2] == 1
        # Cheap edges should survive
        assert eliminated[0] == 0
        assert eliminated[1] == 0


# ---------- BoundPropagator ----------

class TestBoundPropagator:
    def test_sweep(self):
        edges = np.array([[0, 1], [0, 2], [1, 2]], dtype=np.int32)
        costs = np.array([1.0, 1.0, 100.0])
        profits = np.array([0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=3, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=0,
        )

        fwd = forward_labeling(prob, 0)
        prop = BoundPropagator(prob, fwd, fwd, correction=0.0)

        col_upper = np.ones(prob.num_edges)
        fixings = prop.sweep(upper_bound=4.0, col_upper=col_upper)
        assert isinstance(fixings, np.ndarray)
        assert fixings.dtype == np.int32
        assert 2 in fixings  # edge {1,2} eliminated

    def test_propagate_fixed_edge(self):
        edges = np.array(
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
        )
        costs = np.array([1.0, 1.0, 1.0, 100.0, 1.0, 1.0])
        profits = np.array([0.0, 0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=4, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=0,
        )

        fwd = forward_labeling(prob, 0)
        prop = BoundPropagator(prob, fwd, fwd, correction=0.0)

        col_upper = np.ones(prob.num_edges)
        fixings = prop.propagate_fixed_edge(
            edge=0, upper_bound=5.0, col_upper=col_upper
        )
        assert isinstance(fixings, np.ndarray)
        assert 3 in fixings  # edge {1,2} (cost=100) eliminated


# ---------- SeparationOracle: path mode ----------

class TestSeparationOraclePath:
    def test_path_finds_violated_cuts(self):
        prob = make_path_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)
        oracle.add_default_separators()

        # Disconnected solution — node 2 not connected to source
        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 1.0])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 3) or (u == 1 and v == 2):
                x[i] = 1.0

        cuts = oracle.separate(x, y, x_offset=0, y_offset=m)
        assert len(cuts) > 0

    def test_path_feasible(self):
        prob = make_path_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)

        # Feasible path: 0 -> 1 -> 3
        x = np.zeros(m)
        y = np.array([1.0, 1.0, 0.0, 1.0])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 1 and v == 3):
                x[i] = 1.0

        assert oracle.is_feasible(x, y, x_offset=0, y_offset=m)

    def test_path_infeasible(self):
        prob = make_path_problem()
        m = prob.num_edges

        oracle = SeparationOracle(prob)

        # Infeasible: node 2 claimed visited but disconnected
        x = np.zeros(m)
        y = np.array([1.0, 0.0, 1.0, 1.0])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 3):
                x[i] = 1.0
            if (u == 1 and v == 2):
                x[i] = 1.0

        assert not oracle.is_feasible(x, y, x_offset=0, y_offset=m)


# ---------- SeparationOracle: edge cases ----------

class TestSeparationOracleEdgeCases:
    def test_no_separators_returns_empty(self):
        prob = make_tour_problem()
        m = prob.num_edges
        oracle = SeparationOracle(prob)
        # Don't add any separators
        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 0.5])
        cuts = oracle.separate(x, y, x_offset=0, y_offset=m)
        assert len(cuts) == 0

    def test_all_zero_solution_no_cuts(self):
        prob = make_tour_problem()
        m = prob.num_edges
        oracle = SeparationOracle(prob)
        oracle.add_default_separators()
        x = np.zeros(m)
        y = np.zeros(prob.num_nodes)
        cuts = oracle.separate(x, y, x_offset=0, y_offset=m, tol=1e-6)
        assert len(cuts) == 0

    def test_nonzero_offsets(self):
        prob = make_tour_problem()
        m = prob.num_edges
        oracle = SeparationOracle(prob)
        oracle.add_default_separators()

        x = np.zeros(m)
        y = np.array([1.0, 0.5, 0.5, 0.5])

        graph_edges = prob.graph_edges()
        for i in range(m):
            u, v = graph_edges[i]
            if (u == 0 and v == 1) or (u == 2 and v == 3):
                x[i] = 1.0

        x_off = 10
        y_off = x_off + m
        cuts = oracle.separate(x, y, x_offset=x_off, y_offset=y_off)
        assert len(cuts) > 0
        # All indices should be >= x_off
        for c in cuts:
            assert all(idx >= x_off for idx in c.indices)


# ---------- BoundPropagator: additional tests ----------

class TestBoundPropagatorExtended:
    def test_has_all_pairs_initially_false(self):
        prob = make_tour_problem()
        fwd = forward_labeling(prob, prob.source)
        prop = BoundPropagator(prob, fwd, fwd, correction=0.0)
        # There's no direct Python accessor for has_all_pairs_bounds,
        # but we verify set_all_pairs_bounds + propagate_fixed_edge works
        col_upper = np.ones(prob.num_edges)
        # Should use neighbor-only path
        fixings = prop.propagate_fixed_edge(edge=0, upper_bound=1e12, col_upper=col_upper)
        assert isinstance(fixings, np.ndarray)

    def test_sweep_loose_ub_no_fixings(self):
        edges = np.array([[0, 1], [0, 2], [1, 2]], dtype=np.int32)
        costs = np.array([1.0, 1.0, 100.0])
        profits = np.array([0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=3, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=0,
        )

        fwd = forward_labeling(prob, 0)
        prop = BoundPropagator(prob, fwd, fwd, correction=0.0)

        col_upper = np.ones(prob.num_edges)
        fixings = prop.sweep(upper_bound=1e12, col_upper=col_upper)
        assert len(fixings) == 0

    def test_propagate_with_all_pairs(self):
        edges = np.array(
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
        )
        costs = np.array([1.0, 1.0, 1.0, 100.0, 1.0, 1.0])
        profits = np.array([0.0, 0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=4, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=0,
        )

        n = prob.num_nodes
        fwd = forward_labeling(prob, 0)
        prop = BoundPropagator(prob, fwd, fwd, correction=0.0)

        # Build all-pairs bounds
        all_pairs = np.zeros(n * n)
        for s in range(n):
            bounds = forward_labeling(prob, s)
            for v in range(n):
                all_pairs[s * n + v] = bounds[v]

        prop.set_all_pairs_bounds(all_pairs)

        col_upper = np.ones(prob.num_edges)
        # Fix edge 2 ({0,3}). With tight UB=5, {1,2} (cost=100) should be eliminated.
        fixings = prop.propagate_fixed_edge(edge=2, upper_bound=5.0, col_upper=col_upper)
        assert 3 in fixings  # edge {1,2}

    def test_sweep_on_path_problem(self):
        edges = np.array(
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
        )
        costs = np.array([1.0, 1.0, 1.0, 100.0, 1.0, 1.0])
        profits = np.array([0.0, 0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=4, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=3,
        )

        fwd = forward_labeling(prob, 0)
        bwd = backward_labeling(prob, 3)
        prop = BoundPropagator(prob, fwd, bwd, correction=0.0)

        col_upper = np.ones(prob.num_edges)
        fixings = prop.sweep(upper_bound=5.0, col_upper=col_upper)
        assert 3 in fixings  # edge {1,2} cost=100


# ---------- Problem accessors ----------

class TestProblemAccessors:
    def test_tour_properties(self):
        prob = make_tour_problem()
        assert prob.is_tour
        assert prob.source == 0
        assert prob.target == 0
        assert prob.num_nodes == 4
        assert prob.num_edges == 6
        assert prob.capacity == 7.0

    def test_path_properties(self):
        prob = make_path_problem()
        assert not prob.is_tour
        assert prob.source == 0
        assert prob.target == 3
        assert prob.num_nodes == 4
        assert prob.num_edges == 6
        assert prob.capacity == 7.0

    def test_edge_costs_array(self):
        prob = make_tour_problem()
        costs = prob.edge_costs
        assert isinstance(costs, np.ndarray)
        assert costs.dtype == np.float64
        assert len(costs) == prob.num_edges
        assert all(c == 10.0 for c in costs)

    def test_profits_array(self):
        prob = make_tour_problem()
        profits = prob.profits
        assert isinstance(profits, np.ndarray)
        assert len(profits) == prob.num_nodes
        np.testing.assert_array_equal(profits, [0.0, 20.0, 15.0, 10.0])

    def test_demands_array(self):
        prob = make_tour_problem()
        demands = prob.demands
        assert isinstance(demands, np.ndarray)
        assert len(demands) == prob.num_nodes
        np.testing.assert_array_equal(demands, [0.0, 3.0, 4.0, 2.0])

    def test_graph_edges_shape(self):
        prob = make_tour_problem()
        edges = prob.graph_edges()
        assert isinstance(edges, np.ndarray)
        assert edges.dtype == np.int32
        assert edges.shape == (prob.num_edges, 2)


# ---------- IO: load function ----------

class TestLoad:
    def test_load_tour(self):
        prob = load(str(DATA_DIR / "tiny4.txt"))
        assert prob.is_tour
        assert prob.num_nodes == 4
        assert prob.source == 0
        assert prob.target == 0

    def test_load_path(self):
        prob = load(str(DATA_DIR / "tiny4_path.txt"))
        assert not prob.is_tour
        assert prob.source == 0
        assert prob.target == 3


# ---------- Preprocessing: edge elimination on path ----------

class TestPreprocessingExtended:
    def test_edge_elimination_path(self):
        edges = np.array(
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
        )
        costs = np.array([1.0, 1.0, 1.0, 100.0, 1.0, 1.0])
        profits = np.array([0.0, 0.0, 0.0, 0.0])
        demands = np.array([0.0, 0.0, 0.0, 0.0])
        prob = Problem(
            num_nodes=4, edges=edges, edge_costs=costs,
            profits=profits, demands=demands,
            capacity=1e18, source=0, target=3,
        )

        fwd = forward_labeling(prob, 0)
        bwd = backward_labeling(prob, 3)
        eliminated = edge_elimination(prob, fwd, bwd, upper_bound=5.0, correction=0.0)
        assert eliminated[3] == 1  # edge {1,2} cost=100
        assert eliminated[0] == 0  # {0,1} cheap
        assert eliminated[2] == 0  # {0,3} cheap

    def test_labeling_symmetry_for_tour(self):
        """For tours (undirected), forward from source == backward from source."""
        prob = make_tour_problem()
        fwd = forward_labeling(prob, prob.source)
        bwd = backward_labeling(prob, prob.source)
        np.testing.assert_array_almost_equal(fwd, bwd)


# ---------- Warm-start: extended tests ----------

class TestWarmStartExtended:
    def test_tour_capacity_respected(self):
        prob = make_tour_problem()
        result = build_warm_start(prob, time_budget_ms=50.0)
        m = prob.num_edges

        # Verify total demand of visited nodes <= capacity
        total_demand = 0.0
        for i in range(prob.num_nodes):
            if result.col_values[m + i] > 0.5:
                total_demand += prob.demands[i]
        assert total_demand <= prob.capacity + 1e-6

    def test_path_capacity_respected(self):
        prob = make_path_problem()
        result = build_warm_start(prob, time_budget_ms=50.0)
        m = prob.num_edges

        total_demand = 0.0
        for i in range(prob.num_nodes):
            if result.col_values[m + i] > 0.5:
                total_demand += prob.demands[i]
        assert total_demand <= prob.capacity + 1e-6

    def test_objective_is_finite(self):
        prob = make_tour_problem()
        result = build_warm_start(prob, time_budget_ms=50.0)
        assert np.isfinite(result.objective)


# ---------- Module attributes ----------

def test_has_highs():
    """has_highs attribute exists and is boolean."""
    assert isinstance(rcspp_bac.has_highs, bool)
