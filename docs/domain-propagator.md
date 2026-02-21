# Domain Propagator

File: `src/model/highs_bridge.cpp` (install_propagator method)

## Overview

A custom domain propagator injected into HiGHS via the `HighsUserPropagator` patch.
Runs inside `HighsSearch.cpp` after HiGHS's built-in reduced-cost fixing
(`HighsRedcostFixing`). Uses ESPPRC-style labeling bounds to fix edge variables
to zero when no improving tour can use that edge.

## Labeling Bounds

Computed during preprocessing (`src/preprocess/edge_elimination.h`) and passed
to the propagator at setup time.

**Forward labeling** `labeling_from(prob, root)`:
- ESPPRC-style label-correcting algorithm from a given root node
- Label state: (net_cost, demand, predecessor)
- Net cost = sum(edge_costs) - sum(profits) along path
- Seed: f[root] = -profit(root), demand = demand(root)
- 2-cycle elimination: no immediate return to predecessor
- Capacity check: accumulated demand <= Q
- Label dominance: (cost1, demand1) dominates (cost2, demand2) iff
  cost1 <= cost2 AND demand1 <= demand2
- Max 50 labels per node to bound complexity
- Returns f[v] = minimum net cost to reach v from root

Convenience wrappers: `forward_labeling(prob, source)`, `forward_labeling(prob)`,
`backward_labeling(prob, target)`.

**Correction term**: When source = target (tour), the depot profit is subtracted
in both forward and backward bounds. Add `correction = profit(source)` to
compensate. For s-t paths (source != target), correction = 0.

## Trigger A: Upper Bound Sweep

When the MIP upper bound improves (`ub < last_ub - 1e-9`):

1. Reset all processed-edge flags (tighter UB may enable new fixings).
2. For each unfixed edge e = (u, v):
   - Compute lower bound: `lb = min(f[u] + cost(e) + b[v], f[v] + cost(e) + b[u]) + correction`
   - If `lb > ub + 1e-6`: fix x_e = 0 via `domain.changeBound(kUpper, e, 0.0)`
3. For each unfixed node i: if all incident edges are fixed to 0, fix y_i = 0.

This is the same logic as static edge elimination in preprocessing, but applied
dynamically as the upper bound tightens during branch-and-bound.

## Trigger B: Chained Bounds

When an edge x_(a,i) is fixed to 1 (by branching or propagation):

### Default: Neighbor-only scan

Scan edges adjacent to both endpoints of the fixed edge:

**Neighbors of i** (forward direction through fixed edge):
- `cost_a_to_i = f[a] + cost(a,i) - profit(i)`
- For each neighbor edge (i, j): `lb = cost_a_to_i + cost(i,j) + b[j] + correction`
- If `lb > ub + 1e-6`: fix x_(i,j) = 0

**Neighbors of a** (backward direction through fixed edge):
- `cost_via_i_return = cost(a,i) + b[i] - profit(a)`
- For each neighbor edge (k, a): `lb = f[k] + cost(k,a) + cost_via_i_return + correction`
- If `lb > ub + 1e-6`: fix x_(k,a) = 0

### All-pairs variant (--all_pairs_propagation true)

When all-pairs labeling bounds are available, scans **all** unfixed edges (not
just neighbors). For each candidate edge (u, v), tests 8 orientations combining:
- Fixed edge direction: a->i or i->a
- Candidate edge direction: u->v or v->u
- Candidate edge placement: before or after fixed edge in the tour

Each orientation computes: `lb = d[start->...->endpoint1] + cost + d[endpoint2->...->finish] + correction`.

Disabled by default — benchmarks show negligible gain over neighbor-only scan.

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

## CLI Options

| Flag | Default | Description |
|---|---|---|
| `--all_pairs_propagation true` | false | Use all-pairs labeling for stronger Trigger B |

## Preprocessing Pipeline

The propagator reuses labeling bounds computed in `Model::solve()` during the
parallel preprocessing phase:

```
tbb::task_group
├─ forward_labeling(problem, source)    → fwd_bounds
├─ backward_labeling(problem, target)   → bwd_bounds  (same as fwd for tour)
└─ build_warm_start(problem, budget)    → warm_start solution
```

These bounds serve double duty:
1. **Static edge elimination** in `build_formulation()` — fix variables before solve
2. **Dynamic propagation** during solve — fix variables as UB tightens and edges get fixed
