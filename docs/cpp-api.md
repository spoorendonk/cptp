# C++ API

## CMake Integration

Link against `rcspp_model` (which pulls in `rcspp_core` and `rcspp_sep` transitively):

```cmake
# In your CMakeLists.txt
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE rcspp_model)
```

## Solving a Tour (closed loop)

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

## Solving an s-t Path (open)

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

## Loading from File

```cpp
#include "core/io.h"
#include "model/model.h"

auto prob = rcspp::io::load("instance.txt");  // reads source/target from file

rcspp::Model model;
model.set_problem(std::move(prob));
auto result = model.solve({{"time_limit", "120"}});
```

## SolveResult Fields

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
