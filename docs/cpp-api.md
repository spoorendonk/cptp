# C++ API

## CMake Integration

`rcspp-bac` exports one public library target:

```cmake
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE rcspp)
```

HiGHS is mandatory and linked transitively.

## Public Headers

- `core/problem.h`
- `core/io.h`
- `core/solution.h`
- `model/model.h`

## `Model`

```cpp
#include "model/model.h"

rcspp::Model model;
model.set_graph(
    4,
    {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}},
    {10.0, 8.0, 12.0, 6.0, 7.0, 5.0});

model.set_depot(0);            // Tour mode
// model.set_source(0);        // Path mode
// model.set_target(3);

model.set_profits({0.0, 20.0, 15.0, 10.0});
model.add_capacity_resource({0.0, 3.0, 4.0, 2.0}, 7.0);

auto result = model.solve({{"time_limit", "60"}, {"output_flag", "false"}});
```

### Key methods

- `set_problem(Problem)`
- `set_graph(num_nodes, edges, edge_costs)`
- `set_depot(depot)` for tours (`source == target`)
- `set_source(source)` + `set_target(target)` for s-t paths
- `set_profits(profits)`
- `add_capacity_resource(demands, limit)`
- `solve(options)`

## `io::load`

```cpp
#include "core/io.h"
#include "model/model.h"

auto prob = rcspp::io::load("instance.txt");
rcspp::Model model;
model.set_problem(std::move(prob));
auto result = model.solve();
```

## `SolveResult`

- `status`
- `objective`
- `bound`
- `gap`
- `time_seconds`
- `nodes`
- `tour`
- `tour_arcs`
- `separator_stats`
- `total_cuts`
- `separation_rounds`

## Build/Test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/rcspp_tests
```
