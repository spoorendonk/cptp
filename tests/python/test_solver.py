"""Test the Python CPTP solver bindings."""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cptp
from cptp import Model, Problem, SolveResult, Status, SeparatorStats, load, solve

DATA_DIR = Path(__file__).parent.parent / "data"
BIN_PATH = Path(__file__).resolve().parents[2] / "build" / "cptp-solve"
OPTS = [("time_limit", "30"), ("output_flag", "false")]


# ---------- Model.set_graph() API ----------

def test_model_basic():
    """Tour with set_graph: 3 nodes, visits 2 profitable ones."""
    model = Model()
    edges = np.array([[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]], dtype=np.int32)
    costs = np.array([5.0, 5.0, 3.0, 3.0, 4.0, 4.0])

    model.set_graph(3, edges, costs)
    model.set_depot(0)
    model.set_profits(np.array([0.0, 10.0, 8.0]))
    model.add_capacity_resource(np.array([0.0, 2.0, 3.0]), 5.0)

    result = model.solve(OPTS)
    assert result.has_solution()
    assert result.is_optimal()
    assert result.status == Status.Optimal
    assert result.objective <= 0.0


def test_model_negative_costs():
    """Tour with negative edge costs."""
    model = Model()
    edges = np.array([[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]], dtype=np.int32)
    costs = np.array([-2.0, 5.0, 3.0, 3.0, 4.0, -1.0])

    model.set_graph(3, edges, costs)
    model.set_depot(0)
    model.set_profits(np.array([0.0, 5.0, 5.0]))
    model.add_capacity_resource(np.array([0.0, 1.0, 1.0]), 10.0)

    result = model.solve(OPTS)
    assert result.has_solution()
    assert result.objective < 0.0


def test_model_st_path():
    """s-t path mode via set_source/set_target."""
    model = Model()
    edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)
    costs = np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0])

    model.set_graph(4, edges, costs)
    model.set_source(0)
    model.set_target(3)
    model.set_profits(np.array([0.0, 20.0, 15.0, 0.0]))
    model.add_capacity_resource(np.array([0.0, 3.0, 4.0, 0.0]), 7.0)

    result = model.solve(OPTS)
    assert result.has_solution()
    assert result.tour[0] == 0
    assert result.tour[-1] == 3


# ---------- load() + set_problem() ----------

def test_load_numeric_tour():
    """Load numeric tour instance and solve via set_problem."""
    prob = load(str(DATA_DIR / "tiny4.txt"))
    assert prob.num_nodes == 4
    assert prob.num_edges == 6
    assert prob.is_tour
    assert prob.source == prob.target

    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)
    assert result.is_optimal()
    assert result.objective == pytest.approx(-11.0)


def test_load_numeric_path():
    """Load numeric s-t path instance and solve."""
    prob = load(str(DATA_DIR / "tiny4_path.txt"))
    assert not prob.is_tour
    assert prob.source != prob.target

    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)
    assert result.has_solution()
    assert result.tour[0] == prob.source
    assert result.tour[-1] == prob.target


def test_load_tsplib():
    """Load TSPLIB .vrp instance."""
    prob = load(str(DATA_DIR / "tiny3.vrp"))
    assert prob.num_nodes > 0
    assert prob.is_tour

    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)
    assert result.has_solution()


# ---------- Problem class ----------

def test_problem_numpy_accessors():
    """Problem vector accessors return numpy arrays (zero-copy views)."""
    prob = load(str(DATA_DIR / "tiny4.txt"))

    assert isinstance(prob.edge_costs, np.ndarray)
    assert prob.edge_costs.dtype == np.float64
    assert len(prob.edge_costs) == prob.num_edges

    assert isinstance(prob.profits, np.ndarray)
    assert prob.profits.dtype == np.float64
    assert len(prob.profits) == prob.num_nodes

    assert isinstance(prob.demands, np.ndarray)
    assert prob.demands.dtype == np.float64
    assert len(prob.demands) == prob.num_nodes

    edges = prob.graph_edges()
    assert isinstance(edges, np.ndarray)
    assert edges.dtype == np.int32
    assert edges.shape == (prob.num_edges, 2)


def test_problem_constructor():
    """Build a Problem from Python and round-trip through the solver."""
    edges = np.array([[0, 1], [0, 2], [1, 2]], dtype=np.int32)
    costs = np.array([5.0, 3.0, 4.0])
    profits = np.array([0.0, 10.0, 8.0])
    demands = np.array([0.0, 2.0, 3.0])

    prob = Problem(
        num_nodes=3, edges=edges, edge_costs=costs,
        profits=profits, demands=demands,
        capacity=5.0, source=0, target=0, name="test",
    )
    assert prob.num_nodes == 3
    assert prob.num_edges == 3
    assert prob.is_tour
    assert prob.name == "test"
    np.testing.assert_array_equal(prob.edge_costs, costs)
    np.testing.assert_array_equal(prob.profits, profits)

    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)
    assert result.has_solution()


# ---------- SolveResult numpy + stats ----------

def test_result_numpy_tour():
    """SolveResult.tour is a numpy int32 array."""
    prob = load(str(DATA_DIR / "tiny4.txt"))
    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)

    assert isinstance(result.tour, np.ndarray)
    assert result.tour.dtype == np.int32
    assert len(result.tour) > 0

    assert isinstance(result.tour_arcs, np.ndarray)
    assert result.tour_arcs.dtype == np.int32


def test_result_stats():
    """SolveResult exposes cut statistics."""
    prob = load(str(DATA_DIR / "tiny4.txt"))
    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)

    assert isinstance(result.total_cuts, int)
    assert isinstance(result.separation_rounds, int)
    assert isinstance(result.nodes, int)
    assert result.time_seconds >= 0.0
    assert result.gap >= 0.0
    assert isinstance(result.separator_stats, dict)
    for name, stats in result.separator_stats.items():
        assert isinstance(name, str)
        assert isinstance(stats, SeparatorStats)
        assert stats.cuts_added >= 0
        assert stats.rounds_called >= 0
        assert stats.time_seconds >= 0.0


def test_status_enum():
    """Status enum values are accessible."""
    assert Status.Optimal != Status.Error
    assert Status.Feasible != Status.Infeasible
    assert Status.TimeLimit != Status.Unbounded


# ---------- solve() convenience wrapper ----------

def test_solve_wrapper_tour():
    """High-level solve() function for a tour."""
    edges = np.array([[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]], dtype=np.int32)
    result = solve(
        num_nodes=3,
        edges=edges,
        edge_costs=np.array([5.0, 5.0, 3.0, 3.0, 4.0, 4.0]),
        profits=np.array([0.0, 10.0, 8.0]),
        demands=np.array([0.0, 2.0, 3.0]),
        capacity=5.0,
        depot=0,
        time_limit=30.0,
        verbose=False,
    )
    assert result.has_solution()
    assert isinstance(result.tour, np.ndarray)


def test_solve_wrapper_path():
    """High-level solve() with source != target."""
    edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)
    result = solve(
        num_nodes=4,
        edges=edges,
        edge_costs=np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0]),
        profits=np.array([0.0, 20.0, 15.0, 0.0]),
        demands=np.array([0.0, 3.0, 4.0, 0.0]),
        capacity=7.0,
        source=0,
        target=3,
        time_limit=30.0,
    )
    assert result.has_solution()
    assert result.tour[0] == 0
    assert result.tour[-1] == 3


# ---------- CLI (__main__) ----------

def test_cli_tour():
    """python -m cptp on a tour instance."""
    result = subprocess.run(
        [sys.executable, "-m", "cptp", str(DATA_DIR / "tiny4.txt")],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Tour:" in result.stdout


def test_cli_path_override():
    """python -m cptp with --source/--target override."""
    result = subprocess.run(
        [sys.executable, "-m", "cptp", str(DATA_DIR / "tiny4.txt"),
         "--source", "0", "--target", "3"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Path:" in result.stdout


def test_binary_output_has_stage_headers():
    """cptp-solve prints the expected stage headers."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Running HiGHS " in result.stdout
    assert "Lower bounds calculation:" in result.stdout
    assert "fwd_reachable=" in result.stdout
    assert "bwd_reachable=" in result.stdout
    assert "Startup Stage" not in result.stdout
    assert "Construction heuristic:" in result.stdout
    assert "Preprocess:" in result.stdout
    assert "Local search:" in result.stdout
    assert "Preprocess restart:" in result.stdout


def test_binary_output_all_pairs_block_before_local_search():
    """all_pairs_bounds prints a dedicated block before local search."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--all_pairs_bounds", "true",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Preprocess:" in result.stdout
    assert "Lower bounds all pairs calculation:" in result.stdout
    assert "all_pairs_reachable=" in result.stdout
    assert "Local search:" in result.stdout


def test_binary_output_propagator_summary_format():
    """Propagator summary uses sweep+chain fixings and no UB-improvements field."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Propagator:" in result.stdout
    assert "sweep=" in result.stdout
    assert "UB improvements" not in result.stdout


def test_binary_legacy_aliases_rejected():
    """Removed legacy aliases are not silently accepted."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    legacy_flags = [
        "--all_pairs_propagation", "--heuristic_callback",
        "--heuristic_strategy", "--heuristic_budget_ms",
        "--submip_separation", "--disable_heuristics",
        "--max_cuts_per_separator",
        "--min_violation_sec", "--min_violation_rci",
        "--min_violation_multistar", "--min_violation_comb",
        "--min_violation_rglm", "--min_violation_spi",
    ]
    for flag in legacy_flags:
        result = subprocess.run(
            [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
             "--time_limit", "10", "--threads", "1", "--random_seed", "0",
             "--output_flag", "true", flag, "1"],
            capture_output=True, text=True, timeout=60,
        )
        assert result.returncode == 0, f"{flag} caused non-zero exit"
        assert "HiGHS rejected option" in result.stdout, (
            f"{flag} was silently accepted instead of being rejected by HiGHS"
        )


def test_binary_settings_nondefault_params():
    """Settings line reports rc_fixing, cut families, and branching sub-params."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "42",
         "--rc_fixing", "off", "--enable_rglm", "true",
         "--enable_sec", "false", "--max_cuts_rci", "5",
         "--tree_violation_factor", "3",
         "--branch_hyper", "sb", "--branch_hyper_sb_max_depth", "3",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    settings_line = ""
    for line in result.stdout.splitlines():
        if "Settings" in line:
            settings_line = line
            break
    assert settings_line, "No Settings line found"
    assert "rc_fixing=off" in settings_line
    assert "enable_rglm=on" in settings_line
    assert "enable_sec=off" in settings_line
    assert "max_cuts_rci=5" in settings_line
    assert "tree_violation_factor=3.000000" in settings_line
    assert "branch_hyper=sb" in settings_line
    assert "branch_hyper_sb_max_depth=3" in settings_line


def test_binary_settings_separation_tol_and_tree_factor():
    """Settings line reports separation_tol and tree_violation_factor only."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--separation_tol", "0.07",
         "--tree_violation_factor", "10",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    settings_line = ""
    for line in result.stdout.splitlines():
        if "Settings" in line:
            settings_line = line
            break
    assert settings_line, "No Settings line found"
    assert "separation_tol=0.070000" in settings_line
    assert "tree_violation_factor=10.000000" in settings_line
    assert "min_violation_" not in settings_line


def test_binary_tree_violation_factor_negative_clamps_to_zero():
    """Negative tree_violation_factor is clamped to 0 and logged."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--tree_violation_factor", "-2",
         "--output_flag", "true"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "tree_violation_factor=-2 is invalid; clamping to 0" in result.stdout
    settings_line = ""
    for line in result.stdout.splitlines():
        if "Settings" in line:
            settings_line = line
            break
    assert settings_line, "No Settings line found"
    assert "tree_violation_factor=0.000000" in settings_line


def test_binary_output_flag_false_is_silent():
    """cptp-solve with output_flag=false prints no informational output."""
    if not BIN_PATH.exists():
        pytest.skip("cptp-solve binary not built")

    result = subprocess.run(
        [str(BIN_PATH), str(DATA_DIR / "tiny4.txt"),
         "--time_limit", "10", "--threads", "1", "--random_seed", "0",
         "--output_flag", "false"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    # output_flag=false suppresses HiGHS + solver internal logging
    # but CLI still prints Instance header and result summary
    assert "Running HiGHS" not in result.stdout
    assert "Startup Stage" not in result.stdout
