# rcspp-bac

Branch-and-cut solver for a variant of the **Resource Constrained Shortest Path Problem (RCSPP)**.

## Problem Variant

The solver targets an RCSPP with:

- **Single capacity resource**: each vertex has a demand; the total demand on the path must not exceed vehicle capacity Q
- **Vertex profits**: each vertex has a profit collected when visited; the objective minimizes `travel_cost - collected_profit`
- **Optional vertices**: not all vertices need to be visited -- the solver selects the optimal subset
- **Negative-cost cycles**: the graph may contain negative cycles (profits can exceed edge costs), requiring elementary path constraints
- **Tours or s-t paths**: supports closed tours (depot-to-depot) and open s-t paths (source-to-target) with automatic detection

This problem is also known as the **Capacitated Profitable Tour Problem (CPTP)**, introduced by [Jepsen, Petersen, Spoorendonk & Pisinger (2014)](https://doi.org/10.1016/S1572-5286(14)00036-X). The s-t path variant extends it to open paths where source and target differ.

## Features

- **MIP formulation** with binary edge and node variables (undirected)
- **Dynamic cut separation**: subtour elimination (SEC), rounded capacity inequalities (RCI), multistar/GLM inequalities, rounded GLM (RGLM)
- **Gomory-Hu tree** (Gusfield's algorithm) shared across separators for efficient min-cut computation
- **Domain propagator**: labeling-based edge fixing during branch-and-bound
- **Preprocessing**: demand-reachability filtering and ESPPRC-style edge elimination
- **Warm-start heuristic**: parallel greedy construction + local search (2-opt, or-opt, node drop/add) via TBB
- **MIP backend**: HiGHS with custom separator, feasibility check, and propagator callbacks
- **Python bindings** via nanobind (optional)

## Build

Requires GCC 14+ (C++23), CMake 3.25+, and TBB:

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

For Python bindings:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCPTP_BUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

## CLI Usage

```bash
./build/cptp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

Accepts TSPLIB (`.vrp`, `.sppcc`) and PathWyse (`.txt`) instance formats. All options beyond `--source`/`--target` are forwarded to HiGHS (e.g., `--time_limit`, `--threads`, `--output_flag`).

When `source != target`, the solver uses an open s-t path formulation (degree 1 at source/target, degree 2 at intermediates). When `source == target` (default), the standard tour formulation is used.

### Examples

```bash
# Tour (closed loop from depot)
./build/cptp-solve bench/instances/spprclib/B-n45-k6-54.sppcc --time_limit 120

# s-t path (source/target read from file)
./build/cptp-solve tests/data/tiny4_path.txt

# Override source/target via CLI (turns a tour instance into a path)
./build/cptp-solve tests/data/tiny4.txt --source 0 --target 3

# Suppress HiGHS log output
./build/cptp-solve tests/data/tiny4.txt --output_flag false
```

### Instance formats

**PathWyse `.txt`** (used by tests):

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

## C++ API

Link against `cptp_model` (which pulls in `cptp_core` and `cptp_sep` transitively):

```cmake
# In your CMakeLists.txt
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE cptp_model)
```

### Solving a tour (closed loop)

```cpp
#include "model/model.h"

cptp::Model model;

// Undirected graph: 4 nodes, 6 edges
std::vector<cptp::Edge> edges = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};
std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
model.set_graph(4, edges, costs);

model.set_depot(0);  // closed tour from node 0
model.set_profits({0.0, 20.0, 15.0, 10.0});
model.add_capacity_resource({0.0, 3.0, 4.0, 2.0}, /*capacity=*/7.0);

auto result = model.solve({{"time_limit", "60"}});
if (result.has_solution()) {
    // result.tour = ordered node sequence, e.g. [0, 1, 3, 0]
    // result.objective = travel_cost - collected_profit
}
```

### Solving an s-t path (open)

```cpp
#include "model/model.h"

cptp::Model model;

std::vector<cptp::Edge> edges = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};
std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
model.set_graph(4, edges, costs);

// source != target → path mode (degree 1 at endpoints, degree 2 at intermediates)
model.set_source(0);
model.set_target(3);

model.set_profits({0.0, 20.0, 15.0, 10.0});
model.add_capacity_resource({0.0, 3.0, 4.0, 2.0}, /*capacity=*/7.0);

auto result = model.solve({{"time_limit", "60"}, {"output_flag", "false"}});
if (result.has_solution()) {
    // result.tour = [source, ..., target], e.g. [0, 1, 3]
    // result.tour.front() == 0, result.tour.back() == 3
}
```

### Loading from file

```cpp
#include "core/io.h"
#include "model/model.h"

auto prob = cptp::io::load("instance.txt");  // reads source/target from file

cptp::Model model;
model.set_problem(std::move(prob));
auto result = model.solve({{"time_limit", "120"}});
```

### SolveResult fields

| Field | Type | Description |
|---|---|---|
| `status` | `Status` | `Optimal`, `Feasible`, `TimeLimit`, `Infeasible`, `Error` |
| `objective` | `double` | `travel_cost - collected_profit` |
| `bound` | `double` | Best dual bound |
| `gap` | `double` | Relative MIP gap (0.0 = proven optimal) |
| `tour` | `vector<int32_t>` | Ordered node sequence (path or tour) |
| `time_seconds` | `double` | Wall-clock solve time |
| `nodes` | `int64_t` | Branch-and-bound nodes explored |

## Tests

```bash
./build/cptp_tests
```

## Project Structure

```
src/core/        Problem definition, IO (TSPLIB & PathWyse), Dinitz max-flow, Gomory-Hu tree
src/sep/         Cut separators (SEC, RCI, Multistar, RGLM, Comb)
src/model/       HiGHS MIP integration (separators, propagator, callbacks)
src/util/        Utilities (Timer)
src/heuristic/   Warm-start construction + local search
src/preprocess/  Demand-reachability and edge elimination
src/cli/         Command-line solver
python/          nanobind Python bindings
tests/           Catch2 unit and regression tests
docs/            Algorithm documentation
```

## Documentation

- [Algorithms and techniques](docs/algorithms.md) -- formulation, solver pipeline, preprocessing, references
- [Cut separation](docs/separation.md) -- SEC, RCI, Multistar/GLM, RGLM, Comb, cut management
- [Domain propagator](docs/domain-propagator.md) -- labeling-based edge fixing during B&C
- [Warm-start heuristic](docs/warm-start-heuristic.md) -- construction, local search, parallelism

## Dependencies

- [HiGHS](https://highs.dev) v1.10.0 -- MIP solver (fetched via CMake)
- [Catch2](https://github.com/catchorg/Catch2) v3.7.1 -- testing (fetched via CMake)
- [TBB](https://github.com/oneapi-src/oneTBB) 2021.11 -- parallel heuristic

## Benchmarks

Tested against [SPPRCLIB](https://or.rwth-aachen.de/research/spprclib) instances and verified against [PathWyse](https://doi.org/10.1080/10556788.2023.2296978).

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). [A branch-and-cut algorithm for the capacitated profitable tour problem](https://doi.org/10.1016/S1572-5286(14)00036-X). *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Salani, M., Basso, S., & Giuffrida, V. (2024). [PathWyse: a flexible, open-source library for the resource constrained shortest path problem](https://doi.org/10.1080/10556788.2023.2296978). *Optimization Methods and Software*, 39(2).
