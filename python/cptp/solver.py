"""High-level convenience wrapper around the C++ CPTP solver."""

from collections.abc import Sequence

import numpy as np

from cptp._cptp import Model as _Model
from cptp._cptp import SolveResult


def solve(
    num_nodes: int,
    edges: np.ndarray,
    edge_costs: np.ndarray,
    profits: np.ndarray,
    demands: np.ndarray,
    capacity: float | Sequence[float] | np.ndarray,
    depot: int = 0,
    source: int | None = None,
    target: int | None = None,
    time_limit: float = 600.0,
    num_threads: int | None = None,
    verbose: bool = False,
) -> SolveResult:
    """Solve a CPTP instance.

    Args:
        num_nodes: Number of nodes in the graph.
        edges: (m, 2) array of (tail, head) pairs.
        edge_costs: (m,) array of edge costs (can be negative).
        profits: (n,) array of node profits.
        demands: (n,) array of node demands, or an (m, n) array with one row
            per capacity constraint.
        capacity: Vehicle capacity, or m capacities when demands is (m, n).
            A solution satisfies every row: sum_i demands[k, i] y_i <= capacity[k].
        depot: Depot node index (default 0). Sets source = target = depot (tour).
        source: Source node for s-t path. Overrides depot if set.
        target: Target node for s-t path. Overrides depot if set.
        time_limit: Time limit in seconds.
        num_threads: Number of threads.
        verbose: Print solver output.

    Returns:
        SolveResult with tour, objective, gap, etc.
    """
    model = _Model()
    model.set_graph(
        num_nodes,
        np.ascontiguousarray(edges, dtype=np.int32),
        np.ascontiguousarray(edge_costs, dtype=np.float64),
    )

    if source is not None or target is not None:
        model.set_source(source if source is not None else depot)
        model.set_target(target if target is not None else depot)
    else:
        model.set_depot(depot)

    model.set_profits(np.ascontiguousarray(profits, dtype=np.float64))
    demand_rows = np.asarray(demands, dtype=np.float64)
    if demand_rows.ndim == 1:
        if not isinstance(capacity, (int, float, np.integer, np.floating)):
            raise ValueError("capacity must be a number when demands is one-dimensional")
        model.add_capacity_resource(np.ascontiguousarray(demand_rows), float(capacity))
    else:
        capacities = np.atleast_1d(np.asarray(capacity, dtype=np.float64))
        if demand_rows.ndim != 2 or capacities.shape != (demand_rows.shape[0],):
            raise ValueError("demands must be (m, n) with one capacity per row")
        for row, limit in zip(demand_rows, capacities, strict=True):
            model.add_capacity_resource(np.ascontiguousarray(row), float(limit))

    options = [
        ("time_limit", str(time_limit)),
        ("output_flag", "true" if verbose else "false"),
    ]
    if num_threads is not None:
        options.append(("threads", str(num_threads)))
    return model.solve(options)
