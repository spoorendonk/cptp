# C++ API

## CMake Integration

### Full solver (with HiGHS)

Link against `rcspp_model` (which pulls in `rcspp_algorithms` transitively):

```cmake
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE rcspp_model)
```

### Algorithms only (no HiGHS)

Link against `rcspp_algorithms` for separators, propagation, heuristic, and preprocessing:

```cmake
add_subdirectory(path/to/rcspp-bac)
set(RCSPP_BUILD_HIGHS OFF)
target_link_libraries(my_app PRIVATE rcspp_algorithms)
```

Only depends on TBB (`apt install libtbb-dev`).

## Standalone Algorithm API

These components are solver-independent and can be used from any MIP solver's cut callback.

### SeparationOracle

Header: `sep/separation_oracle.h`

Bundles support-graph construction, Gomory-Hu tree computation, and parallel separator execution into a single `separate()` call.

```cpp
#include "sep/separation_oracle.h"

rcspp::Problem prob = rcspp::io::load("instance.txt");

rcspp::sep::SeparationOracle oracle(prob);
oracle.add_default_separators();  // SEC + RCI + Multistar + Comb
oracle.set_max_cuts_per_separator(3);

// Inside your solver's cut callback:
// x_values = edge variable LP values (size = num_edges)
// y_values = node variable LP values (size = num_nodes)
// x_offset, y_offset = column offsets in your solver's LP
auto cuts = oracle.separate(x_values, y_values, x_offset, y_offset);

for (const auto& cut : cuts) {
    // cut.indices — column indices (int32_t)
    // cut.values  — coefficients (double)
    // cut.rhs     — right-hand side (double, <= form)
    // cut.violation — how much the current LP violates this cut
    // Add to your solver's cut pool...
}
```

Individual separators can be added selectively:

```cpp
rcspp::sep::SeparationOracle oracle(prob);
oracle.add_separator(std::make_unique<rcspp::sep::SECSeparator>());
oracle.add_separator(std::make_unique<rcspp::sep::RCISeparator>());
// Skip Multistar and Comb if not needed
```

#### Feasibility check (lazy constraints)

Check whether an integer solution satisfies all subtour elimination constraints:

```cpp
if (!oracle.is_feasible(x_int, y_int, x_offset, y_offset)) {
    // Reject incumbent — subtour violation found
}
```

This uses a tight tolerance (1e-6) and only checks SECs. Useful as a lazy-constraint callback on integer-feasible solutions from solver heuristics.

### BoundPropagator

Header: `preprocess/bound_propagator.h`

Domain propagation using labeling bounds. Two triggers:

- **Trigger A (sweep)**: When the upper bound improves, scan all edges for fixings.
- **Trigger B (chain)**: When an edge is fixed to 1, infer fixings for neighboring edges.

```cpp
#include "preprocess/bound_propagator.h"
#include "preprocess/edge_elimination.h"

auto fwd = rcspp::preprocess::forward_labeling(prob, prob.source());
auto bwd = rcspp::preprocess::backward_labeling(prob, prob.target());
double correction = prob.is_tour() ? prob.profits()[prob.source()] : 0.0;

rcspp::preprocess::BoundPropagator prop(prob, fwd, bwd, correction);

// Trigger A: sweep all edges when UB improves
// col_upper[e] = current upper bound of edge variable e
auto fixed_edges = prop.sweep(upper_bound, col_upper);

// Trigger B: chain fixings when edge e is branched to 1
auto chain_fixings = prop.propagate_fixed_edge(e, upper_bound, col_upper);

// After sweep, find nodes with all edges fixed to 0
auto fixed_nodes = prop.sweep_nodes(col_upper, /*y_offset=*/num_edges);
```

For stronger Trigger B propagation with all-pairs bounds:

```cpp
auto dist = rcspp::preprocess::all_pairs_labeling(prob);
prop.set_all_pairs_bounds(std::move(dist));
// Now propagate_fixed_edge scans all edges (not just neighbors)
```

### Warm-Start Heuristic

Header: `heuristic/warm_start.h`

Parallel greedy construction + local search (2-opt, or-opt, node drop/add):

```cpp
#include "heuristic/warm_start.h"

auto warm = rcspp::heuristic::build_warm_start(prob, /*budget_ms=*/500.0);
// warm.col_values — solution vector (size = num_edges + num_nodes)
// warm.objective  — travel_cost - collected_profit

// Feed to your solver as an initial solution
```

### Preprocessing

Header: `preprocess/edge_elimination.h`

Static edge elimination before solving:

```cpp
#include "preprocess/edge_elimination.h"

auto fwd = rcspp::preprocess::forward_labeling(prob, prob.source());
auto bwd = rcspp::preprocess::backward_labeling(prob, prob.target());
double correction = prob.is_tour() ? prob.profits()[prob.source()] : 0.0;

// Returns vector<bool>: eliminated[e] = true if edge e can be removed
auto eliminated = rcspp::preprocess::edge_elimination(
    prob, fwd, bwd, upper_bound, correction);
```

### Testing the Algorithm API

The standalone algorithm API is tested by `rcspp_algo_tests` (57 Catch2 tests, no HiGHS required):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_HIGHS=OFF
cmake --build build -j$(nproc)
./build/rcspp_algo_tests
```

Tests cover:
- **SeparationOracle** (13 tests): cut generation on violated/feasible solutions, path mode, non-zero offsets, cut struct validation, max-cuts limits, cumulative separator addition
- **Individual separators** (11 tests): SEC (tour + path, target-in-S handling), RCI, Multistar, Comb, RGLM, separator name accessors
- **BoundPropagator** (11 tests): Trigger A sweep, Trigger B chain fixings, all-pairs bounds, path mode, sweep_nodes, edge skipping, accessor correctness
- **Preprocessing** (8 tests): demand reachability (tour + path), edge elimination, labeling bounds
- **Warm-start heuristic** (6 tests): tour/path construction, degree/capacity consistency
- **Core** (5 tests): Dinitz max-flow, Gomory-Hu tree, Problem accessors, IO parsing

Run individual test tags:

```bash
./build/rcspp_algo_tests [oracle]       # SeparationOracle tests
./build/rcspp_algo_tests [propagator]   # BoundPropagator tests
./build/rcspp_algo_tests [sec]          # SEC separator tests
./build/rcspp_algo_tests [path]         # s-t path mode tests
./build/rcspp_algo_tests [heuristic]    # Warm-start tests
```

## Solver API (requires HiGHS)

### Solving a Tour (closed loop)

```cpp
#include "model/model.h"

rcspp::Model model;

// Undirected graph: 4 nodes, 6 edges
std::vector<rcspp::Edge> edges = {
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

### Solving an s-t Path (open)

```cpp
#include "model/model.h"

rcspp::Model model;

std::vector<rcspp::Edge> edges = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};
std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
model.set_graph(4, edges, costs);

// source != target -> path mode (degree 1 at endpoints, degree 2 at intermediates)
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

### Loading from File

```cpp
#include "core/io.h"
#include "model/model.h"

auto prob = rcspp::io::load("instance.txt");  // reads source/target from file

rcspp::Model model;
model.set_problem(std::move(prob));
auto result = model.solve({{"time_limit", "120"}});
```

### SolveResult Fields

| Field | Type | Description |
|---|---|---|
| `status` | `Status` | `Optimal`, `Feasible`, `TimeLimit`, `Infeasible`, `Error` |
| `objective` | `double` | `travel_cost - collected_profit` |
| `bound` | `double` | Best dual bound |
| `gap` | `double` | Relative MIP gap (0.0 = proven optimal) |
| `tour` | `vector<int32_t>` | Ordered node sequence (path or tour) |
| `time_seconds` | `double` | Wall-clock solve time |
| `nodes` | `int64_t` | Branch-and-bound nodes explored |
| `total_cuts` | `int` | Total user cuts added |
| `separation_rounds` | `int` | Number of separation rounds |
| `separator_stats` | `map<string, SeparatorStats>` | Per-separator statistics |

Methods: `has_solution()`, `is_optimal()`.
