# Python API

## Installation

```bash
pip install rcspp-bac
```

Or build from source:

```bash
pip install .
```

## Quick Start with `solve()`

### Tour (closed loop)

```python
import numpy as np
from rcspp_bac import solve

edges = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]], dtype=np.int32)
costs = np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0])
profits = np.array([0.0, 20.0, 15.0, 10.0])
demands = np.array([0.0, 3.0, 4.0, 2.0])

result = solve(
    num_nodes=4,
    edges=edges,
    edge_costs=costs,
    profits=profits,
    demands=demands,
    capacity=7.0,
    depot=0,
    time_limit=60.0,
)

if result.has_solution():
    print(result.tour)       # e.g. [0, 1, 3, 0]
    print(result.objective)  # travel_cost - collected_profit
```

### s-t Path (open)

```python
result = solve(
    num_nodes=4,
    edges=edges,
    edge_costs=costs,
    profits=profits,
    demands=demands,
    capacity=7.0,
    source=0,
    target=3,
    time_limit=60.0,
)

if result.has_solution():
    print(result.tour)  # e.g. [0, 1, 3] -- source first, target last
```

### `solve()` Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `num_nodes` | `int` | required | Number of nodes |
| `edges` | `ndarray[int32]` | required | `(m, 2)` array of `(tail, head)` pairs |
| `edge_costs` | `ndarray[float64]` | required | `(m,)` edge costs (can be negative) |
| `profits` | `ndarray[float64]` | required | `(n,)` node profits |
| `demands` | `ndarray[float64]` | required | `(n,)` node demands |
| `capacity` | `float` | required | Vehicle capacity |
| `depot` | `int` | `0` | Depot node (tour mode). Ignored if `source`/`target` set. |
| `source` | `int \| None` | `None` | Source node for s-t path |
| `target` | `int \| None` | `None` | Target node for s-t path |
| `time_limit` | `float` | `600.0` | Time limit in seconds |
| `num_threads` | `int \| None` | `None` | Number of threads |
| `verbose` | `bool` | `False` | Print HiGHS solver output |

## `Model` Class

For finer control, use `Model` directly:

```python
from rcspp_bac import Model
import numpy as np

model = Model()
model.set_graph(
    4,
    np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]], dtype=np.int32),
    np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0]),
)

model.set_depot(0)  # tour mode
# Or for s-t path:
# model.set_source(0)
# model.set_target(3)

model.set_profits(np.array([0.0, 20.0, 15.0, 10.0]))
model.add_capacity_resource(np.array([0.0, 3.0, 4.0, 2.0]), 7.0)

result = model.solve([("time_limit", "60"), ("output_flag", "false")])
```

## Loading from File

```python
from rcspp_bac import Model, load

problem = load("instance.txt")  # auto-detects format

model = Model()
model.set_problem(problem)
result = model.solve([("time_limit", "120")])
```

The `Problem` object exposes instance data as numpy arrays:

```python
problem.num_nodes      # int
problem.num_edges      # int
problem.source         # int
problem.target         # int
problem.capacity       # float
problem.is_tour        # bool
problem.edge_costs     # ndarray[float64]
problem.profits        # ndarray[float64]
problem.demands        # ndarray[float64]
problem.graph_edges()  # ndarray[int32] shape (m, 2)
```

## CLI Usage

```bash
python -m rcspp_bac instance.txt --time_limit 120 --verbose
python -m rcspp_bac instance.sppcc --source 0 --target 5
```

## SolveResult Fields

| Field | Type | Description |
|---|---|---|
| `status` | `Status` | `Optimal`, `Feasible`, `TimeLimit`, `Infeasible`, `Error` |
| `objective` | `float` | `travel_cost - collected_profit` |
| `bound` | `float` | Best dual bound |
| `gap` | `float` | Relative MIP gap (0.0 = proven optimal) |
| `tour` | `ndarray[int32]` | Ordered node sequence (path or tour) |
| `time_seconds` | `float` | Wall-clock solve time |
| `nodes` | `int` | Branch-and-bound nodes explored |
| `total_cuts` | `int` | Total user cuts added |
| `separation_rounds` | `int` | Number of separation rounds |
| `separator_stats` | `dict[str, SeparatorStats]` | Per-separator statistics |

Methods: `has_solution()`, `is_optimal()`.
