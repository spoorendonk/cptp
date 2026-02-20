"""High-level convenience wrapper around the C++ CPTP solver."""

import numpy as np
from cptp._cptp import Model as _Model, SolveResult


def solve(
    num_nodes: int,
    edges: np.ndarray,
    edge_costs: np.ndarray,
    profits: np.ndarray,
    demands: np.ndarray,
    capacity: float,
    depot: int = 0,
    time_limit: float = 600.0,
    num_threads: int = 1,
    verbose: bool = False,
) -> SolveResult:
    """Solve a CPTP instance.

    Args:
        num_nodes: Number of nodes in the graph.
        edges: (m, 2) array of (tail, head) pairs.
        edge_costs: (m,) array of edge costs (can be negative).
        profits: (n,) array of node profits.
        demands: (n,) array of node demands.
        capacity: Vehicle capacity.
        depot: Depot node index (default 0).
        time_limit: Time limit in seconds.
        num_threads: Number of threads.
        verbose: Print solver output.

    Returns:
        SolveResult with tour, objective, gap, etc.
    """
    model = _Model()
    model.set_graph(num_nodes, edges.astype(np.int32), edge_costs.astype(np.float64))
    model.set_depot(depot)
    model.set_profits(profits.astype(np.float64))
    model.add_capacity_resource(demands.astype(np.float64), capacity)
    return model.solve(time_limit=time_limit, num_threads=num_threads, verbose=verbose)
