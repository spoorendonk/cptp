# rcspp-bac

Branch-and-cut solver for the **Resource Constrained Shortest Path Problem (RCSPP)**.

## Scope

This repository now exposes a **solver-only public interface**:

- C++: `rcspp::Model` + core types (`Problem`, `SolveResult`, `Status`)
- Python: `Model`, `solve`, `Problem`, `SolveResult`, `Status`, `load`
- CLI: `rcspp-solve`

Cut separators, preprocessing, propagator, and heuristic modules remain internal implementation details.

## Build

Requirements: GCC 14+, CMake 3.25+, TBB.

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

HiGHS is required and fetched automatically by CMake.

## C++ Usage

Link against the single library target `rcspp`:

```cmake
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE rcspp)
```

Minimal example:

```cpp
#include "model/model.h"

rcspp::Model model;
model.set_graph(4,
    {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}},
    {10.0, 8.0, 12.0, 6.0, 7.0, 5.0});
model.set_depot(0);
model.set_profits({0.0, 20.0, 15.0, 10.0});
model.add_capacity_resource({0.0, 3.0, 4.0, 2.0}, 7.0);

auto result = model.solve({{"time_limit", "60"}});
```

## Python Usage

```bash
pip install .
```

```python
import numpy as np
from rcspp_bac import solve

result = solve(
    num_nodes=4,
    edges=np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]], dtype=np.int32),
    edge_costs=np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0]),
    profits=np.array([0.0, 20.0, 15.0, 10.0]),
    demands=np.array([0.0, 3.0, 4.0, 2.0]),
    capacity=7.0,
    depot=0,
)
```

## CLI

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

## Tests

```bash
./build/rcspp_tests
pytest tests/python/test_solver.py
RCSPP_RUN_LOCAL_INSTALL_TEST=1 pytest tests/python/test_local_install.py
```

## Docs

- [Getting Started](docs/getting-started.md)
- [C++ API](docs/cpp-api.md)
- [Python API](docs/python-api.md)
- [Instance Formats](docs/instance-formats.md)
- [Algorithms Overview](docs/algorithms.md)
- [Cut Separation](docs/separation.md)
- [Domain Propagator](docs/domain-propagator.md)
- [Primal Heuristic](docs/primal-heuristic.md)
- [Preprocessing](docs/preprocessing.md)
