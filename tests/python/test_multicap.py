"""m capacity rows from Python (#90): 2-D demands, capacities, solve()."""

import numpy as np
import pytest
from cptp import Model, Problem, solve

OPTS = [("time_limit", "30"), ("output_flag", "false")]

# Depot 0 and customers 1..3 on a complete graph with unit costs. Row 0 alone
# admits all three customers (-29); row 1 cuts that to {1, 3} (-20); row 2
# forbids 1 with 3, leaving {2, 3} (-18). Same instance as test_multicap.cpp.
EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)
COSTS = np.ones(6)
PROFITS = np.array([0.0, 12.0, 10.0, 11.0])
ROWS = np.array(
    [[0.0, 1.0, 1.0, 1.0], [0.0, 2.0, 2.0, 1.0], [0.0, 1.0, 0.0, 1.0]],
)
CAPS = np.array([3.0, 3.0, 1.0])


def test_problem_two_dimensional_demands() -> None:
    prob = Problem(
        num_nodes=4,
        edges=EDGES,
        edge_costs=COSTS,
        profits=PROFITS,
        demands=ROWS,
        capacities=CAPS,
        source=0,
        target=0,
    )
    assert prob.num_capacities == 3
    np.testing.assert_array_equal(prob.capacities, CAPS)
    np.testing.assert_array_equal(prob.demand_rows, ROWS)
    np.testing.assert_array_equal(prob.demands, ROWS[0])
    assert prob.capacity == CAPS[0]
    # demand_rows is a copy the caller owns.
    rows = prob.demand_rows
    rows[1, 1] = 99.0
    assert prob.demand_rows[1, 1] == 2.0

    model = Model()
    model.set_problem(prob)
    result = model.solve(OPTS)
    assert result.is_optimal()
    assert result.objective == pytest.approx(-18.0)


def test_problem_rejects_mismatched_rows() -> None:
    with pytest.raises(ValueError, match="one capacity per demand row"):
        Problem(
            num_nodes=4,
            edges=EDGES,
            edge_costs=COSTS,
            profits=PROFITS,
            demands=ROWS,
            capacities=CAPS[:2],
            source=0,
            target=0,
        )


@pytest.mark.parametrize(("rows", "objective"), [(1, -29.0), (2, -20.0), (3, -18.0)])
def test_solve_with_m_rows(rows: int, objective: float) -> None:
    result = solve(
        num_nodes=4,
        edges=EDGES,
        edge_costs=COSTS,
        profits=PROFITS,
        demands=ROWS[:rows],
        capacity=list(CAPS[:rows]),
        time_limit=30.0,
    )
    assert result.is_optimal()
    assert result.objective == pytest.approx(objective)
    visited = set(int(v) for v in result.tour)
    for k in range(rows):
        assert sum(ROWS[k, v] for v in visited) <= CAPS[k] + 1e-9


def test_solve_one_dimensional_demands_unchanged() -> None:
    result = solve(
        num_nodes=4,
        edges=EDGES,
        edge_costs=COSTS,
        profits=PROFITS,
        demands=ROWS[1],
        capacity=float(CAPS[1]),
        time_limit=30.0,
    )
    assert result.is_optimal()
    assert result.objective == pytest.approx(-20.0)
