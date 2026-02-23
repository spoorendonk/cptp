# Python API

## Installation

```bash
pip install rcspp-bac
```

Or build from source:

```bash
pip install .
```

The package has two tiers:

- **Algorithm API** (always available) — `SeparationOracle`, `BoundPropagator`, `build_warm_start`, labeling, edge elimination
- **Solver API** (requires HiGHS) — `Model`, `solve()`, CLI

Check at runtime:

```python
import rcspp_bac
if rcspp_bac.has_highs:
    from rcspp_bac import Model, solve
```

## Algorithm API (no HiGHS required)

These are solver-independent — use them from any MIP solver's cut callback (Gurobi, CPLEX, SCIP, etc.).

### SeparationOracle

Cut separation for LP relaxation solutions:

```python
import numpy as np
from rcspp_bac import Problem, SeparationOracle

prob = Problem(
    num_nodes=4,
    edges=np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]], dtype=np.int32),
    edge_costs=np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0]),
    profits=np.array([0.0, 20.0, 15.0, 10.0]),
    demands=np.array([0.0, 3.0, 4.0, 2.0]),
    capacity=7.0, source=0, target=0,
)

oracle = SeparationOracle(prob)
oracle.add_default_separators()  # SEC + RCI + Multistar + Comb
oracle.set_max_cuts_per_separator(3)

# Inside your solver's cut callback:
# x_values = edge variable LP values, y_values = node variable LP values
cuts = oracle.separate(x_values, y_values, x_offset=0, y_offset=num_edges)

for cut in cuts:
    print(cut.indices)    # ndarray[int32] — column indices
    print(cut.values)     # ndarray[float64] — coefficients
    print(cut.rhs)        # float — right-hand side (<= form)
    print(cut.violation)  # float — how much the LP violates this cut
    print(cut.size)       # int — number of nonzeros
    # Add to your solver's cut pool...
```

Individual separators can be added selectively:

```python
oracle = SeparationOracle(prob)
oracle.add_sec()        # Subtour Elimination Constraints
oracle.add_rci()        # Rounded Capacity Inequalities
oracle.add_multistar()  # Multistar / GLM
oracle.add_comb()       # Comb inequalities
oracle.add_rglm()       # Rounded GLM
```

#### Feasibility check (lazy constraints)

```python
# Check if an integer solution has subtour violations
feasible = oracle.is_feasible(x_int, y_int, x_offset=0, y_offset=num_edges)
if not feasible:
    # Reject incumbent — add violated SECs as lazy constraints
    pass
```

### Warm-Start Heuristic

Parallel greedy construction + local search:

```python
from rcspp_bac import build_warm_start

warm = build_warm_start(prob, time_budget_ms=500.0)
warm.col_values  # ndarray[float64] — solution vector (edges + nodes)
warm.objective   # float — travel_cost - collected_profit
```

The solution vector has size `num_edges + num_nodes`:
- `col_values[e] = 1` for edges in the route
- `col_values[num_edges + v] = 1` for visited nodes

### Preprocessing

Labeling-based edge elimination and bound computation:

```python
from rcspp_bac import forward_labeling, backward_labeling, edge_elimination

# Compute labeling bounds
fwd = forward_labeling(prob, prob.source)   # ndarray[float64], size = num_nodes
bwd = backward_labeling(prob, prob.target)  # ndarray[float64], size = num_nodes

# Static edge elimination
correction = 0.0  # or prob.profits[prob.source] for tours
eliminated = edge_elimination(prob, fwd, bwd, upper_bound=100.0, correction=correction)
# eliminated: ndarray[int8], 1 = can be removed, 0 = keep
```

### BoundPropagator

Domain propagation for branch-and-bound integration:

```python
from rcspp_bac import BoundPropagator

prop = BoundPropagator(prob, fwd, bwd, correction=0.0)

# Trigger A: sweep all edges when UB improves
col_upper = np.ones(prob.num_edges)
fixed_edges = prop.sweep(upper_bound=100.0, col_upper=col_upper)
# fixed_edges: ndarray[int32] — edge indices that can be fixed to 0

# Trigger B: chain fixings when edge e is branched to 1
chain = prop.propagate_fixed_edge(edge=0, upper_bound=100.0, col_upper=col_upper)
# chain: ndarray[int32] — additional edges that can be fixed to 0
```

### Problem

Build a problem from Python or load from file:

```python
from rcspp_bac import Problem, load

# From file (auto-detects TSPLIB or numeric format)
prob = load("instance.txt")

# From Python arrays
prob = Problem(
    num_nodes=4,
    edges=np.array([[0,1], [0,2], [1,2]], dtype=np.int32),
    edge_costs=np.array([5.0, 3.0, 4.0]),
    profits=np.array([0.0, 10.0, 8.0]),
    demands=np.array([0.0, 2.0, 3.0]),
    capacity=5.0,
    source=0,     # source == target => tour mode
    target=0,
    name="my_instance",
)
```

Properties (read-only, zero-copy numpy views):

```python
prob.num_nodes      # int
prob.num_edges      # int
prob.source         # int
prob.target         # int
prob.capacity       # float
prob.is_tour        # bool (source == target)
prob.edge_costs     # ndarray[float64]
prob.profits        # ndarray[float64]
prob.demands        # ndarray[float64]
prob.graph_edges()  # ndarray[int32] shape (m, 2)
```

## Solver API (requires HiGHS)

### Quick Start with `solve()`

#### Tour (closed loop)

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

#### s-t Path (open)

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

#### `solve()` Parameters

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

### `Model` Class

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

### Loading from File

```python
from rcspp_bac import Model, load

problem = load("instance.txt")  # auto-detects format

model = Model()
model.set_problem(problem)
result = model.solve([("time_limit", "120")])
```

### CLI Usage

```bash
python -m rcspp_bac instance.txt --time_limit 120 --verbose
python -m rcspp_bac instance.sppcc --source 0 --target 5
```

### SolveResult Fields

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
