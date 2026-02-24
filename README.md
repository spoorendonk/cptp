# rcspp-bac

Branch-and-cut solver for a variant of the **Resource Constrained Shortest Path Problem (RCSPP)**.

> **Disclaimer:** This project was developed entirely through [Claude Code](https://docs.anthropic.com/en/docs/claude-code). 😱

## Problem Variant

The solver targets an RCSPP with:

- **Single capacity resource**: each vertex has a demand; the total demand on the path must not exceed vehicle capacity Q
- **Vertex profits**: each vertex has a profit collected when visited; the objective minimizes `travel_cost - collected_profit`
- **Optional vertices**: not all vertices need to be visited -- the solver selects the optimal subset
- **Negative-cost cycles**: the graph may contain negative cycles (profits can exceed edge costs), requiring elementary path constraints
- **Tours or s-t paths**: supports closed tours (depot-to-depot) and open s-t paths (source-to-target) with automatic detection

This problem generalizes the **Capacitated Profitable Tour Problem (CPTP)** of [Jepsen et al. (2014)](https://doi.org/10.1016/S1572-5286(14)00036-X), which is the special case where source equals target (closed tour).

## Use Cases

The primary application is as a **pricing solver in branch-and-price for CVRP**. In column generation for vehicle routing, the pricing subproblem is an RCSPP where dual variables create vertex profits and negative-cost cycles. Standard labeling algorithms (Pessoa, Sadykov, Uchoa & Vanderbeck, 2020) solve the pricing subproblem via dynamic programming; this solver provides a branch-and-cut alternative following the approach of Jepsen et al. (2014), revisited by Paro (2022) with modern open-source MIP solvers.

## Features

- **Pluggable algorithm library** (`rcspp_algorithms`) — separators, propagation, heuristic, and preprocessing, usable with any MIP solver
- **Dynamic cut separation**: subtour elimination (SEC), rounded capacity inequalities (RCI), multistar/GLM inequalities, rounded GLM (RGLM), comb
- **Gomory-Hu tree** (Gusfield's algorithm) shared across separators for efficient min-cut computation
- **Domain propagator**: labeling-based edge fixing during branch-and-bound
- **Preprocessing**: demand-reachability filtering and labeling-based edge elimination
- **Warm-start heuristic**: parallel greedy construction + local search (2-opt, or-opt, node drop/add) via TBB
- **Batteries-included HiGHS integration** — full branch-and-cut solver out of the box
- **Python bindings** via nanobind (optional), with algorithm API available without HiGHS

## Architecture

The project is split into two layers:

1. **`rcspp_algorithms`** — solver-independent library containing all cutting-plane separators, bound propagation, warm-start heuristic, and preprocessing. Only depends on TBB. Can be used standalone from any MIP solver's cut callback.

2. **`rcspp_model`** (optional) — HiGHS integration that wires the algorithms into HiGHS's branch-and-cut framework. Three callback interfaces are patched into HiGHS (see `third_party/highs_patch/`):
   - **User separator** (`HighsUserSeparator`) — cut separation at fractional and integer-feasible nodes
   - **Feasibility checker** — lazy constraint validation rejecting subtour-violating incumbents
   - **Domain propagator** (`HighsUserPropagator`) — labeling-based variable fixing during B&B

## Build

Requires GCC 14+ (C++23), CMake 3.25+, and TBB:

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

To build only the algorithm library (no HiGHS dependency):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_HIGHS=OFF
cmake --build build -j$(nproc)
```

For Python bindings:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

## CLI Usage

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

Accepts TSPLIB (`.vrp`, `.sppcc`) and numeric (`.txt`) instance formats. All options beyond `--source`/`--target`/`--branch_hyper` are forwarded to HiGHS (e.g., `--time_limit`, `--threads`, `--output_flag`).

#### Hyperplane branching

The `--branch_hyper` option enables dynamic constraint branching, which adds/removes LP constraint rows during the branch-and-bound search. Modes:

| Mode | Description |
|------|-------------|
| `off` (default) | Standard variable branching only |
| `pairs` | Ryan-Foster pairs: branch on `y_i + y_j` for nearest-neighbor node pairs |
| `clusters` | Cluster demand: branch on `sum(d_i * y_i)` for small node clusters |
| `demand` | Global demand: branch on total demand `sum(d_i * y_i)` |
| `cardinality` | Cardinality: branch on `sum(y_i)` (number of visited nodes) |
| `all` | All of the above combined |

When `source != target`, the solver uses an open s-t path formulation (degree 1 at source/target, degree 2 at intermediates). When `source == target` (default), the standard tour formulation is used.

### Examples

```bash
# Tour (closed loop from depot)
./build/rcspp-solve benchmarks/instances/spprclib/B-n45-k6-54.sppcc --time_limit 120

# s-t path (source/target read from file)
./build/rcspp-solve tests/data/tiny4_path.txt

# Override source/target via CLI (turns a tour instance into a path)
./build/rcspp-solve tests/data/tiny4.txt --source 0 --target 3

# Suppress HiGHS log output
./build/rcspp-solve tests/data/tiny4.txt --output_flag false

# Hyperplane branching (Ryan-Foster pairs)
./build/rcspp-solve bench/instances/spprclib/B-n45-k6-54.sppcc --branch_hyper pairs
```

### Instance formats

**Numeric `.txt`** (used by tests):

```
4 6               # num_nodes num_edges
0 1 10            # u v cost  (repeated for each edge)
0 2 8
0 3 12
1 2 6
1 3 7
2 3 5
7                 # capacity Q
0 3               # source target (optional; omit or use "0 0" for tour)
```

Node profits and demands are read from the same format (see `src/core/io.cpp` for full spec). **TSPLIB `.sppcc`/`.vrp`** files are also supported with the standard `DEPOT_SECTION` (always treated as tour).

## Standalone Algorithm API

The separators, heuristic, and propagator can be used independently of HiGHS — plug them into any MIP solver's cut callback (Gurobi, CPLEX, SCIP, etc.).

### C++ (link against `rcspp_algorithms`)

```cpp
#include "sep/separation_oracle.h"
#include "preprocess/bound_propagator.h"
#include "heuristic/warm_start.h"

// Build or load a Problem
rcspp::Problem prob = rcspp::io::load("instance.txt");

// --- Cut separation ---
rcspp::sep::SeparationOracle oracle(prob);
oracle.add_default_separators();  // SEC + RCI + Multistar + Comb

// Inside your solver's cut callback:
auto cuts = oracle.separate(x_values, y_values, x_offset, y_offset);
for (auto& cut : cuts) {
    // cut.indices, cut.values, cut.rhs — add to your solver
}

// --- Feasibility check (lazy constraints) ---
if (!oracle.is_feasible(x_int, y_int, x_offset, y_offset)) {
    // reject incumbent — subtour violation
}

// --- Warm-start heuristic ---
auto warm = rcspp::heuristic::build_warm_start(prob, /*budget_ms=*/500.0);
// warm.col_values — initial solution, warm.objective — objective value

// --- Bound propagation ---
auto fwd = rcspp::preprocess::forward_labeling(prob, prob.source());
auto bwd = rcspp::preprocess::backward_labeling(prob, prob.target());
double correction = prob.is_tour() ? prob.profits()[prob.source()] : 0.0;

rcspp::preprocess::BoundPropagator prop(prob, fwd, bwd, correction);
auto fixed_edges = prop.sweep(upper_bound, col_upper);       // Trigger A
auto chain = prop.propagate_fixed_edge(e, upper_bound, col_upper);  // Trigger B
```

### Python (no HiGHS required)

```python
import numpy as np
from rcspp_bac import (
    Problem, SeparationOracle, BoundPropagator,
    build_warm_start, forward_labeling, backward_labeling, edge_elimination,
)

prob = Problem(num_nodes=4, edges=edges, edge_costs=costs,
               profits=profits, demands=demands, capacity=7.0)

# Cut separation
oracle = SeparationOracle(prob)
oracle.add_default_separators()
cuts = oracle.separate(x_values, y_values, x_offset=0, y_offset=num_edges)

for cut in cuts:
    # cut.indices (int32), cut.values (float64), cut.rhs — add to solver
    pass

# Warm-start heuristic
warm = build_warm_start(prob, time_budget_ms=500.0)

# Bound propagation
fwd = forward_labeling(prob, prob.source)
bwd = backward_labeling(prob, prob.target)
prop = BoundPropagator(prob, fwd, bwd, correction=0.0)
fixed = prop.sweep(upper_bound=100.0, col_upper=np.ones(prob.num_edges))
```

See [C++ API](docs/cpp-api.md) and [Python API](docs/python-api.md) for full documentation.

## Solver API (HiGHS)

- **Python**: `pip install rcspp-bac` — [usage and examples](docs/python-api.md)
- **C++**: link against `rcspp_model` — [usage and examples](docs/cpp-api.md)

## Tests

```bash
./build/rcspp_algo_tests                    # 57 C++ algorithm tests (no HiGHS needed)
./build/rcspp_tests                         # 26 C++ integration tests (requires HiGHS)
pytest tests/python/test_algorithms.py      # 33 Python algorithm tests (no HiGHS needed)
pytest tests/python/test_solver.py          # Python solver tests (requires HiGHS)
RCSPP_RUN_LOCAL_INSTALL_TEST=1 pytest tests/python/test_local_install.py  # venv + pip install -e . smoke test
```

The algorithm tests (`rcspp_algo_tests`) cover all solver-independent components:
- **SeparationOracle**: cut generation, feasibility checks, offset handling, cut limits, tour and path modes
- **Individual separators**: SEC, RCI, Multistar, Comb, RGLM — both violation detection and correctness
- **BoundPropagator**: Trigger A sweep, Trigger B chain fixings, all-pairs bounds, path mode
- **Preprocessing**: demand reachability, edge elimination, labeling bounds
- **Warm-start heuristic**: tour/path construction, degree/capacity consistency
- **Core**: Dinitz max-flow, Gomory-Hu tree, Problem accessors, IO parsing

Python algorithm tests mirror the C++ coverage for the `rcspp_bac` package bindings.

## Project Structure

```
src/core/        Problem definition, IO (TSPLIB & numeric), Dinitz max-flow, Gomory-Hu tree
src/sep/         Cut separators (SEC, RCI, Multistar, RGLM, Comb) + SeparationOracle
src/preprocess/  BoundPropagator, demand-reachability, edge elimination
src/heuristic/   Warm-start construction + local search
src/model/       HiGHS MIP integration (optional — separators, propagator, callbacks)
src/cli/         Command-line solver (requires HiGHS)
src/util/        Utilities (Logger, Timer)
python/          nanobind Python bindings
tests/           Catch2 unit tests + Python test suites
benchmarks/      Benchmark instances, scripts, and results
docs/            Algorithm documentation
```

### Libraries

| Target | Depends on | Description |
|---|---|---|
| `rcspp_algorithms` | TBB | Solver-independent: separators, propagation, heuristic, preprocessing |
| `rcspp_model` | `rcspp_algorithms` + HiGHS | HiGHS integration (optional, controlled by `RCSPP_BUILD_HIGHS`) |

## Documentation

- [Python API](docs/python-api.md) — installation, algorithm API, `solve()`, `Model`, CLI
- [C++ API](docs/cpp-api.md) — CMake integration, standalone algorithms, solver examples
- [Instance formats](docs/instance-formats.md) — TSPLIB, numeric `.txt`
- [Algorithms and techniques](docs/algorithms.md) — formulation, solver pipeline, references
- [Cut separation](docs/separation.md) — SEC, RCI, Multistar/GLM, RGLM, Comb, SeparationOracle
- [Preprocessing](docs/preprocessing.md) — demand-reachability filtering, labeling-based edge elimination
- [Domain propagator](docs/domain-propagator.md) — BoundPropagator, labeling-based edge fixing
- [Warm-start heuristic](docs/warm-start-heuristic.md) — construction, local search, parallelism
- [Benchmark results](docs/benchmarks.md) — SPPRCLIB and Roberti instance results

## Dependencies

- [TBB](https://github.com/oneapi-src/oneTBB) — parallel heuristic and separators (required)
- [HiGHS](https://highs.dev) — MIP solver (fetched via CMake, optional with `-DRCSPP_BUILD_HIGHS=OFF`)
- [Catch2](https://github.com/catchorg/Catch2) — testing (fetched via CMake)
- [nanobind](https://github.com/wjakob/nanobind) — Python bindings (fetched via CMake, optional)

## Benchmarks

Tested on 76 instances from two standard sets: [SPPRCLIB](https://or.rwth-aachen.de/research/spprclib) (45 instances) and Roberti Set 3 from [Jepsen et al. (2014)](https://doi.org/10.1016/S1572-5286(14)00036-X) (31 instances). **71 of 76 solved to proven optimality** (93%) within a 1-hour time limit. The 5 unsolved are large M-series instances (151--200 nodes).

See [docs/benchmarks.md](docs/benchmarks.md) for full results with UB/LB, gap, timing, and per-separator cut counts.

## Related Projects

- [PathWyse](https://doi.org/10.1080/10556788.2023.2296978) -- labeling-based exact solver for RCSPP (Salani, Basso & Giuffrida 2024)
- [dparo/master-thesis](https://github.com/dparo/master-thesis) -- branch-and-cut pricer for CVRP using CPLEX (Paro 2022)

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). [A branch-and-cut algorithm for the capacitated profitable tour problem](https://doi.org/10.1016/S1572-5286(14)00036-X). *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Pessoa, A., Sadykov, R., Uchoa, E., & Vanderbeck, F. (2020). A generic exact solver for vehicle routing and related problems. *Mathematical Programming*, 183, 483-523.
- Paro, D. (2022). [Exact algorithms for capacitated vehicle routing problems](https://github.com/dparo/master-thesis). Master's thesis, University of Padua.
- Salani, M., Basso, S., & Giuffrida, V. (2024). [PathWyse: a flexible, open-source library for the resource constrained shortest path problem](https://doi.org/10.1080/10556788.2023.2296978). *Optimization Methods and Software*, 39(2).
