"""High-level convenience wrapper around the C++ RCSPP solver."""

import numpy as np
from rcspp_bac._rcspp_bac import Model as _Model, SolveResult


def solve(
    num_nodes: int,
    edges: np.ndarray,
    edge_costs: np.ndarray,
    profits: np.ndarray,
    demands: np.ndarray,
    capacity: float,
    depot: int = 0,
    source: int | None = None,
    target: int | None = None,
    time_limit: float = 600.0,
    num_threads: int = 1,
    verbose: bool = False,
) -> SolveResult:
    """Solve an RCSPP instance.

    Args:
        num_nodes: Number of nodes in the graph.
        edges: (m, 2) array of (tail, head) pairs.
        edge_costs: (m,) array of edge costs (can be negative).
        profits: (n,) array of node profits.
        demands: (n,) array of node demands.
        capacity: Vehicle capacity.
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
    model.set_graph(num_nodes, edges.astype(np.int32), edge_costs.astype(np.float64))

    if source is not None or target is not None:
        model.set_source(source if source is not None else depot)
        model.set_target(target if target is not None else depot)
    else:
        model.set_depot(depot)

    model.set_profits(profits.astype(np.float64))
    model.add_capacity_resource(demands.astype(np.float64), capacity)

    options = [
        ("time_limit", str(time_limit)),
        ("threads", str(num_threads)),
        ("output_flag", "true" if verbose else "false"),
    ]
    return model.solve(options)
