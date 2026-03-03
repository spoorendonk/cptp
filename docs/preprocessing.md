# Preprocessing

Files: `src/preprocess/reachability.h`, `src/preprocess/edge_elimination.h`, `src/preprocess/ng_labeling.h`

## Overview

Preprocessing is staged before the MIP is built:

1. parallel: fast source/target 2-cycle labeling (forward/backward) + fast warm-start
2. edge-elimination trigger computation from Stage-1 bounds and fast UB
3. parallel second warm-start (optional, adaptive), all-pairs 2-cycle
   propagation (optional), and s-t ng/DSSR refinement
4. hand final `(fwd, bwd, UB)` to model build and propagator

An async ng/DSSR worker can continue tightening during branch-and-bound.

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
- `ng_max_size` (default 6)
- `ng_dssr_iters` (default 6)
- `ng_simd` (default true)
- `preproc_adaptive` (default true)
- `preproc_fast_restarts` (default 12)
- `preproc_fast_budget_ms` (default 30)
- `preproc_second_ws_large_n` (default 60)
- `preproc_second_ws_min_elim` (default 0.05)
- `preproc_second_ws_min_elim_large` (default 0.02)
- `preproc_second_ws_budget_ms_min` (default 20)
- `preproc_second_ws_budget_ms_max` (default 400)
- `preproc_second_ws_budget_scale` (default 8)

## Edge Elimination

File: `src/preprocess/edge_elimination.h`

Given forward bounds `f` and backward bounds `b`:

- **Tour**: `lb(u,v) = f[u] + c(u,v) + f[v] + profit(source)`
- **Path**: `lb(u,v) = min(f_s[u] + c(u,v) + f_t[v], f_s[v] + c(u,v) + f_t[u])`

If `lb(u,v) > UB`, edge `(u,v)` is fixed to zero before MIP build.

## Async DSSR During B&B

Files: `src/preprocess/shared_bounds.h`, `src/model/model.cpp`,
`src/model/highs_bridge.cpp`

If `dssr_background_updates=true`, a background worker runs tighter fixed-ng
labeling stages while HiGHS branch-and-bound is running:

- `dssr_background_policy=fixed` runs all planned stages (subject to budget/caps)
- `dssr_background_policy=auto` stops early when there is no observed search
  checkpoint activity or no recent DSSR-stage improvements
- `dssr_background_max_epochs=N` hard-caps the number of async DSSR stages
  (`0` means uncapped)
- `dssr_background_auto_min_epochs=N` sets the minimum stages before `auto`
  early-stop checks (default 4)
- `dssr_background_auto_no_progress_limit=N` sets consecutive non-improving
  stage limit for `auto` (default 6)

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
