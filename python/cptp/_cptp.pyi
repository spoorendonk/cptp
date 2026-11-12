"""CPTP Branch-and-Cut Solver"""

from collections.abc import Sequence
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


class Status(enum.Enum):
    Optimal = 0

    Feasible = 1

    Infeasible = 2

    Unbounded = 3

    TimeLimit = 4

    Error = 5

class SeparatorStats:
    @property
    def cuts_added(self) -> int: ...

    @property
    def rounds_called(self) -> int: ...

    @property
    def time_seconds(self) -> float: ...

class SolveResult:
    @property
    def status(self) -> Status: ...

    @property
    def objective(self) -> float: ...

    @property
    def bound(self) -> float: ...

    @property
    def gap(self) -> float: ...

    @property
    def time_seconds(self) -> float: ...

    @property
    def nodes(self) -> int: ...

    @property
    def tour(self) -> Annotated[NDArray[numpy.int32], dict(shape=(None,), writable=False)]: ...

    @property
    def tour_arcs(self) -> Annotated[NDArray[numpy.int32], dict(shape=(None,), writable=False)]: ...

    @property
    def total_cuts(self) -> int: ...

    @property
    def separation_rounds(self) -> int: ...

    @property
    def separator_stats(self) -> dict[str, SeparatorStats]: ...

    def is_optimal(self) -> bool: ...

    def has_solution(self) -> bool: ...

class Problem:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, num_nodes: int, edges: Annotated[NDArray[numpy.int32], dict(shape=(None, 2))], edge_costs: Annotated[NDArray[numpy.float64], dict(shape=(None,))], profits: Annotated[NDArray[numpy.float64], dict(shape=(None,))], demands: Annotated[NDArray[numpy.float64], dict(shape=(None,))], capacity: float, source: int = 0, target: int = 0, name: str = '') -> None: ...

    @overload
    def __init__(self, num_nodes: int, edges: Annotated[NDArray[numpy.int32], dict(shape=(None, 2))], edge_costs: Annotated[NDArray[numpy.float64], dict(shape=(None,))], profits: Annotated[NDArray[numpy.float64], dict(shape=(None,))], demands: Annotated[NDArray[numpy.float64], dict(shape=(None, None))], capacities: Annotated[NDArray[numpy.float64], dict(shape=(None,))], source: int = 0, target: int = 0, name: str = '') -> None:
        """Build with m capacity rows: demands is (m, n), capacities is (m,)"""

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def num_nodes(self) -> int: ...

    @property
    def num_edges(self) -> int: ...

    @property
    def source(self) -> int: ...

    @property
    def target(self) -> int: ...

    @property
    def capacity(self) -> float:
        """Capacity of row 0 (the only row when num_capacities == 1)"""

    @property
    def num_capacities(self) -> int:
        """Number of capacity rows m; a solution satisfies all of them"""

    @property
    def capacities(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), writable=False)]:
        """(m,) capacities Q_k, a read-only view"""

    @property
    def is_tour(self) -> bool: ...

    @property
    def edge_costs(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), writable=False)]: ...

    @property
    def profits(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), writable=False)]: ...

    @property
    def demands(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), writable=False)]:
        """(n,) demands of row 0, a read-only view"""

    @property
    def demand_rows(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """(m, n) demands d^k_i of every row, a copy"""

    def graph_edges(self) -> Annotated[NDArray[numpy.int32], dict(shape=(None, 2))]:
        """Return (m, 2) numpy array of (tail, head) pairs"""

class Model:
    def __init__(self) -> None: ...

    def set_problem(self, problem: Problem) -> None: ...

    def problem(self) -> Problem: ...

    def set_graph(self, num_nodes: int, edges: Annotated[NDArray[numpy.int32], dict(shape=(None, 2))], edge_costs: Annotated[NDArray[numpy.float64], dict(shape=(None,))]) -> None: ...

    def set_source(self, arg: int, /) -> None: ...

    def set_target(self, arg: int, /) -> None: ...

    def set_depot(self, arg: int, /) -> None: ...

    def set_profits(self, profits: Annotated[NDArray[numpy.float64], dict(shape=(None,))]) -> None: ...

    def add_capacity_resource(self, demands: Annotated[NDArray[numpy.float64], dict(shape=(None,))], limit: float) -> None:
        """
        Add a capacity row sum_i demands[i] * y_i <= limit; every row added is enforced
        """

    def solve(self, options: Sequence[tuple[str, str]] = []) -> SolveResult: ...

has_highs: bool = True

def load(path: str) -> Problem:
    """Load an CPTP instance from file (auto-detect format)"""
