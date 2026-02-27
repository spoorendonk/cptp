# API (C++ and Python)

## Python Install

```bash
pip install rcspp-bac
```

From source:

```bash
pip install .
```

## C++

CMake target:

```cmake
add_subdirectory(path/to/rcspp-bac)
target_link_libraries(my_app PRIVATE rcspp)
```

Minimal usage:

```cpp
#include "model/model.h"

rcspp::Model model;
model.set_graph(4,
    {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}},
    {10.0, 8.0, 12.0, 6.0, 7.0, 5.0});
model.set_depot(0);
model.set_profits({0.0, 20.0, 15.0, 10.0});
model.add_capacity_resource({0.0, 3.0, 4.0, 2.0}, 7.0);

auto result = model.solve({{"time_limit", "60"}, {"output_flag", "false"}});
```

## Python

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
    time_limit=60.0,
)
```

For path mode, pass `source` and `target`.

## Core Types

- `Problem`
- `Model`
- `SolveResult`
- `Status`
- `load`
- `solve` (Python convenience wrapper)
