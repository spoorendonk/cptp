from rcspp_bac._rcspp_bac import (
    # Core types (always available)
    Problem,
    SeparatorStats,
    SolveResult,
    Status,
    load,
    # Solver-independent algorithms
    Cut,
    SeparationOracle,
    WarmStartResult,
    BoundPropagator,
    build_warm_start,
    forward_labeling,
    backward_labeling,
    edge_elimination,
    has_highs,
)

__all__ = [
    # Core types
    "Problem",
    "SeparatorStats",
    "SolveResult",
    "Status",
    "load",
    # Algorithm API
    "Cut",
    "SeparationOracle",
    "WarmStartResult",
    "BoundPropagator",
    "build_warm_start",
    "forward_labeling",
    "backward_labeling",
    "edge_elimination",
    "has_highs",
]

# HiGHS-dependent imports (optional)
if has_highs:
    from rcspp_bac._rcspp_bac import Model
    from rcspp_bac.solver import solve
    __all__ += ["Model", "solve"]
