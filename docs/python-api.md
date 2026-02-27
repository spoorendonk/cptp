# Python API

## Installation

```bash
pip install rcspp-bac
```

Or from source:

```bash
pip install .
```

## Public API

```python
from rcspp_bac import (
    Problem,
    Model,
    solve,
    SolveResult,
    Status,
    load,
    has_highs,
)
```

`has_highs` is always `True` in supported builds.

## Quick Start: `solve()`

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

## `Model`

```python
import numpy as np
from rcspp_bac import Model

model = Model()
model.set_graph(
    4,
    np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]], dtype=np.int32),
    np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0]),
)
model.set_depot(0)
model.set_profits(np.array([0.0, 20.0, 15.0, 10.0]))
model.add_capacity_resource(np.array([0.0, 3.0, 4.0, 2.0]), 7.0)

result = model.solve([("time_limit", "60"), ("output_flag", "false")])
```

## `Problem` and `load()`

```python
from rcspp_bac import load

prob = load("instance.txt")
```

You can also construct `Problem` directly from numpy arrays.

## Tests

```bash
pytest tests/python/test_solver.py
RCSPP_RUN_LOCAL_INSTALL_TEST=1 pytest tests/python/test_local_install.py
```
