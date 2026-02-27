# Domain Propagator

Files: `src/preprocess/bound_propagator.h` (solver-independent), `src/model/highs_bridge.cpp` (HiGHS integration)

## Overview

Labeling-based domain propagation fixes edge variables to zero when no improving
tour/path can use that edge. Available in two forms:

1. **`BoundPropagator`** (solver-independent) — standalone class in `src/preprocess/bound_propagator.h`, usable from any MIP solver's propagation callback.
2. **HiGHS integration** — injected via `HighsUserPropagator` patch, runs inside `HighsSearch.cpp` after HiGHS's built-in reduced-cost fixing.

## BoundPropagator (solver-independent)

File: `src/preprocess/bound_propagator.h`, `src/preprocess/bound_propagator.cpp`

Encapsulates labeling bounds, adjacency structures, and both propagation triggers
into a reusable class. Pre-builds edge costs, profits, and adjacency lists in the
constructor for efficient repeated calls.

```cpp
#include "preprocess/bound_propagator.h"
#include "preprocess/edge_elimination.h"

auto fwd = rcspp::preprocess::forward_labeling(prob, prob.source());
auto bwd = rcspp::preprocess::backward_labeling(prob, prob.target());
double correction = prob.is_tour() ? prob.profits()[prob.source()] : 0.0;

rcspp::preprocess::BoundPropagator prop(prob, fwd, bwd, correction);

// Trigger A: sweep all edges when UB improves
auto fixed = prop.sweep(upper_bound, col_upper);
auto fixed_nodes = prop.sweep_nodes(col_upper, /*y_offset=*/num_edges);

// Trigger B: chain fixings when edge is branched to 1
auto chain = prop.propagate_fixed_edge(edge, upper_bound, col_upper);

// Optional: all-pairs bounds for stronger Trigger B
prop.set_all_pairs_bounds(dist);  // flat n*n array
```

### Methods

| Method | Returns | Description |
|---|---|---|
| `sweep(ub, col_upper)` | `vector<int32_t>` | Trigger A: edge indices fixable to 0 |
| `sweep_nodes(col_upper, y_offset)` | `vector<int32_t>` | Node variable indices fixable to 0 (after sweep) |
| `propagate_fixed_edge(e, ub, col_upper)` | `vector<int32_t>` | Trigger B: edges fixable to 0 given edge e=1 |
| `set_all_pairs_bounds(dist)` | void | Enable stronger all-pairs Trigger B |

## HiGHS Integration

The HiGHS propagator uses the same algorithmic logic as `BoundPropagator` but is
implemented as an inline lambda in `HiGHSBridge::install_propagator()` for
performance (avoids vector allocations in the hot path). It interacts directly
with HiGHS's `HighsDomain` for bound changes.

## Labeling Bounds

Computed during preprocessing via `src/preprocess/ng_labeling.h` and passed to
the propagator at setup time.

**ng-route / DSSR labeling** (`ng::compute_bounds`):
- Label-correcting algorithm with iterative ng-neighborhood growth
- Label state: $(\text{net\_cost},\ \text{demand},\ \text{predecessor},\ \text{ng-visited bitset})$
- Net cost = $\sum \text{edge\_costs} - \sum \text{profits}$ along path
- Capacity check: accumulated $\text{demand} \le Q$
- Dominance uses cost, demand, and visited-subset relation
- SIMD prefilter path (AVX2) for candidate scans with scalar fallback
- Returns forward/backward lower bounds used by both edge elimination and propagation

When `dssr_async=true`, additional bound snapshots are published during solve
through `SharedBoundsStore`, and the propagator callback switches to newer
versions without restarting HiGHS.

**Correction term**: When source = target (tour), the depot profit is subtracted
in both forward and backward bounds. Add $\text{correction} = \text{profit}(\text{source})$ to
compensate. For s-t paths (source $\ne$ target), $\text{correction} = 0$.

## Trigger A: Upper Bound Sweep

When the MIP upper bound improves ($\text{ub} < \text{last\_ub} - 10^{-9}$):

1. Reset all processed-edge flags (tighter UB may enable new fixings).
2. For each unfixed edge $e = (u, v)$:
   - Compute lower bound: $\text{lb} = \min(f[u] + \text{cost}(e) + b[v],\ f[v] + \text{cost}(e) + b[u]) + \text{correction}$
   - If $\text{lb} > \text{ub} + 10^{-6}$: fix $x_e = 0$ via `domain.changeBound(kUpper, e, 0.0)`
3. For each unfixed node $i$: if all incident edges are fixed to 0, fix $y_i = 0$.

This is the same logic as static edge elimination in preprocessing, but applied
dynamically as the upper bound tightens during branch-and-bound.

## Trigger B: Chained Bounds

When an edge $x_{(a,i)}$ is fixed to 1 (by branching or propagation):

### Default: Neighbor-only scan

Scan edges adjacent to both endpoints of the fixed edge:

**Neighbors of i** (forward direction through fixed edge):
- $\text{cost}_{a \to i} = f[a] + \text{cost}(a,i) - \text{profit}(i)$
- For each neighbor edge $(i, j)$: $\text{lb} = \text{cost}_{a \to i} + \text{cost}(i,j) + b[j] + \text{correction}$
- If $\text{lb} > \text{ub} + 10^{-6}$: fix $x_{(i,j)} = 0$

**Neighbors of a** (backward direction through fixed edge):
- $\text{cost\_via\_i\_return} = \text{cost}(a,i) + b[i] - \text{profit}(a)$
- For each neighbor edge $(k, a)$: $\text{lb} = f[k] + \text{cost}(k,a) + \text{cost\_via\_i\_return} + \text{correction}$
- If $\text{lb} > \text{ub} + 10^{-6}$: fix $x_{(k,a)} = 0$

### All-pairs variant (--all_pairs_propagation true)

When all-pairs labeling bounds are available, scans **all** unfixed edges (not
just neighbors). For each candidate edge (u, v), tests 8 orientations combining:
- Fixed edge direction: a->i or i->a
- Candidate edge direction: u->v or v->u
- Candidate edge placement: before or after fixed edge in the tour

Each orientation computes: `lb = d[start->...->endpoint1] + cost + d[endpoint2->...->finish] + correction`.

Disabled by default — benchmarks show negligible gain over neighbor-only scan.

## Testing

The BoundPropagator is covered by 11 C++ tests (`[propagator]` tag) and 5 Python tests:

| Test | Description |
|---|---|
| Sweep fixes expensive edges | Trigger A eliminates cost=100 edge with tight UB |
| Sweep with loose UB | Confirms no fixings with very large UB |
| Sweep skips already-fixed edges | Pre-fixed edges not reported again |
| Sweep on s-t path problem | Path-mode Trigger A with bidirectional bounds |
| sweep_nodes fixes isolated nodes | Nodes with all incident edges fixed to 0 |
| sweep_nodes no fixings when edges open | All edges open → no node fixings |
| propagate_fixed_edge (neighbor-only) | Trigger B eliminates expensive neighbor |
| propagate_fixed_edge (all-pairs) | Trigger B with all-pairs bounds eliminates distant edges |
| propagate_fixed_edge loose UB | Confirms no chain fixings with large UB |
| has_all_pairs_bounds initially false | Verifies default state |
| Accessor correctness | fwd_bounds, bwd_bounds, correction values |

```bash
./build/rcspp_algo_tests [propagator]       # All BoundPropagator tests
./build/rcspp_algo_tests [propagator][path] # Path-mode propagator tests
```

## HiGHS Integration

### Patching

The propagator requires two HiGHS source patches:

- `HighsUserPropagator.h/.cpp`: Callback interface with `setCallback()` /
  `clearCallback()` / `propagate()` static methods. Copied to HiGHS `src/mip/`
  via CMake PATCH_COMMAND.
- `HighsSearch.cpp`: Modified to call `HighsUserPropagator::propagate()` after
  `HighsRedcostFixing`.

### Shared state

The callback lambda captures shared_ptrs for mutable state that persists across
invocations:
- `last_ub`: tracks previous upper bound for Trigger A detection
- `processed_edge`: bitmap of edges already processed by Trigger B
- Statistics counters: propagator_calls, ub_improvements, sweep_fixings, chain_fixings

### Sub-MIP guard

The callback checks `mipsolver.submip` and returns immediately for sub-MIPs.
Sub-MIPs have a different column space — using original model column indices
would crash.

### Statistics

Reported at solve completion:
```
Propagator: 1234 calls, 5 UB improvements, 89 fixings (42 sweep + 47 chain)
```

## Trigger C: Lagrangian Reduced-Cost Fixing (edges to 0)

Uses LP reduced costs as edge weights in capacity-aware labeling to obtain
Lagrangian bounds, then fixes edges whose Lagrangian excess proves they cannot
be in any improving solution.

1. Extract reduced costs from the LP relaxation: `rc_edge[e] = col_dual[e]`,
   `rc_profit[i] = -col_dual[m+i]`.
2. Run forward (and backward for s-t paths) labeling with RC costs.
3. Compute `z_LR` = cheapest resource-feasible tour/path under RC costs.
   - Tour: `z_LR = min_e (rc_fwd[u] + rc_e + rc_fwd[v]) + correction`
   - Path: `z_LR = rc_fwd[target]`
4. For each unfixed edge e=(u,v):
   - `excess(e) = min(rc_fwd[u] + rc_e + rc_bwd[v], rc_fwd[v] + rc_e + rc_bwd[u]) + corr - z_LR`
   - If `z_LP + excess(e) > UB + 1e-6`: fix x_e = 0
5. Fix isolated nodes (all incident edges fixed to 0) to y_i = 0.

This is strictly stronger than Trigger A because it combines LP dual information
with resource-feasibility constraints from the labeling.

## Trigger D: Lagrangian Reduced-Cost Fixing (nodes to 1)

Optional (enabled with `--rc_fixing_to_one true`). For each unfixed node i, runs
labeling with node i forbidden (`profit_i = -1e30`). If the Lagrangian gap from
excluding node i exceeds the UB, the node must be visited in every optimal solution.

- `z_LR_{-i}` = cheapest tour/path avoiding node i
- If `z_LP + (z_LR_{-i} - z_LR) > UB + 1e-6`: fix y_i = 1

Requires N labeling runs per invocation (parallelized via TBB). Disabled by
default as benchmarks show minimal benefit (only 6 total node fixings across 28
SPPRCLIB instances).

## CLI Options

| Flag | Default | Description |
|---|---|---|
| `--rc_fixing <strategy>` | `adaptive` | Lagrangian RC fixing: `off`, `root_only`, `on_ub_improvement`, `periodic`, `adaptive` |
| `--rc_fixing_interval N` | 100 | Interval for `periodic` strategy |
| `--rc_fixing_to_one true` | false | Enable Trigger D (fix nodes to 1) |
| `--all_pairs_propagation true` | false | Use all-pairs labeling for stronger Trigger B |
| `--dssr_async true` | true | Run background DSSR bound tightening during B&B |
| `--ng_initial_size N` | 4 | Initial ng neighborhood size |
| `--ng_max_size N` | 12 | Maximum ng neighborhood size |
| `--ng_dssr_iters N` | 6 | DSSR iterations per labeling run |
| `--ng_label_budget N` | 50 | Max active labels per node in ng labeling |
| `--ng_simd true` | true | Enable AVX2 prefilter path for dominance candidate scans |

## Preprocessing Pipeline

The propagator reuses labeling bounds computed in `Model::solve()` during the
parallel preprocessing phase:

```
tbb::task_group
├─ ng::compute_bounds(problem, source, target, opts)  → fwd_bounds, bwd_bounds
└─ build_initial_solution(problem, budget)            → initial solution
```

These bounds serve double duty:
1. **Static edge elimination** in `build_formulation()` — fix variables before solve
2. **Dynamic propagation** during solve — fix variables as UB tightens and edges get fixed

With async DSSR enabled, a background worker publishes tighter snapshots while
branch-and-bound is running, and the same Trigger A/B logic consumes the
refreshed arrays.

When an async ng/DSSR run finds a full feasible path, it is published as an
external incumbent candidate and injected via the MIP user-solution callback.
If its objective matches/exceeds the current dual proof condition, HiGHS is
interrupted and the solve is reported as optimal when primal and dual bounds
coincide.

`adaptive` RC fixing policy runs on UB-improvement events and self-disables
after repeated low-yield rounds.
