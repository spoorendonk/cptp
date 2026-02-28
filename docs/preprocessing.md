# Preprocessing

Files: `src/preprocess/reachability.h`, `src/preprocess/edge_elimination.h`, `src/preprocess/ng_labeling.h`

## Overview

Two preprocessing phases run in parallel (via TBB) before the MIP is built:

1. demand-reachability node filtering
2. ng-route/DSSR labeling bounds for edge elimination and propagation

The warm-start heuristic runs in parallel with these phases and provides the
initial upper bound used by elimination.

## Demand-Reachability

File: `src/preprocess/reachability.h`

Eliminates nodes that are unreachable within the capacity budget.

- **Tour** (`source == target`): Dijkstra from depot; eliminate nodes where
  `2 * dist(depot, v) - demand(v) > Q`
- **Path** (`source != target`): Dijkstra from source and target; eliminate
  nodes where `dist_s(v) + dist_t(v) - demand(v) > Q`

## ng-route / DSSR Labeling

File: `src/preprocess/ng_labeling.h`

The solver computes forward/backward lower-bound arrays with iterative ng-path
labeling (DSSR style):

- Label state tracks `(cost, demand, predecessor, ng-visited bitset)`
- Dominance checks `(cost, demand, visited-subset)`
- AVX2 prefilter is used when available (`ng_simd=true`)
- ng neighborhoods are grown iteratively from detected cycles until reaching a
  configured cap

Main options (from `Model::solve` options map):

- `ng_initial_size` (default 4)
- `ng_max_size` (default 12)
- `ng_dssr_iters` (default 6)
- `ng_label_budget` (default 50)
- `ng_simd` (default true)

## Edge Elimination

File: `src/preprocess/edge_elimination.h`

Given forward bounds `f` and backward bounds `b`:

- **Tour**: `lb(u,v) = f[u] + c(u,v) + f[v] + profit(source)`
- **Path**: `lb(u,v) = min(f_s[u] + c(u,v) + f_t[v], f_s[v] + c(u,v) + f_t[u])`

If `lb(u,v) > UB`, edge `(u,v)` is fixed to zero before MIP build.

## Async DSSR During B&B

Files: `src/preprocess/shared_bounds.h`, `src/model/model.cpp`,
`src/model/highs_bridge.cpp`

If `dssr_async=true`, a background worker runs tighter fixed-ng labeling stages
while HiGHS branch-and-bound is running:

- publishes improved `(f,b)` snapshots via `SharedBoundsStore`
- snapshots are monotone tightened element-wise
- propagator callback consumes newer versions on the fly

This allows stronger pruning without restarting the solve.

## Testing

Relevant coverage:

- labeling and propagation regression tests in `tests/test_propagator.cpp`
- preprocessing/propagation behavior exercised through `Model::solve()` integration tests
  in `tests/test_model.cpp` and `tests/test_propagator.cpp`
- end-to-end solver validation in Python tests (`tests/python/test_solver.py`)

Run:

```bash
./build/rcspp_tests
pytest tests/python/test_solver.py
```
