# Primal Heuristic

File: `src/heuristic/primal_heuristic.h` (header-only)

## Purpose

Two heuristic components:
1. **Initial solution** (`build_initial_solution`): runs before MIP solve on the complete graph, provides an initial feasible solution to HiGHS via `setSolution()`. Solver-independent — can be used standalone or fed to any MIP solver as a warm start.
2. **LP-guided callback** (`lp_guided_heuristic`): runs during MIP solve via `kCallbackMipUserSolution`, builds reduced graphs from LP relaxation values and runs construction + local search on them

## Algorithm

### Construction (greedy insertion)

Given an ordering of customers, insert each one at the cheapest position in the route:

1. Initialize:
   - **Tour** (source == target): route = [depot, depot]
   - **Path** (source != target): route = [source, target]
2. For each customer c in the given order:
   - Find position p in route minimizing `cost(prev,c) + cost(c,next) - cost(prev,next)`
   - Insert if capacity allows
3. If too few customers inserted, force-insert the cheapest (tours need at least 2 customers for distinct edges)

The insertion order determines which customers get added first (and consume capacity).
Three deterministic orderings are used:
- **Profit/demand ratio** (descending) -- best value per unit capacity
- **Absolute profit** (descending)
- **Distance from source** (ascending) -- short detour customers first

Additional restarts use random permutations.

### Local search

After construction, run improving moves until no improvement found (max 200 iterations):

| Neighborhood | Description |
|---|---|
| **2-opt** | Reverse a segment of the route |
| **Or-opt** | Relocate a chain of 1, 2, or 3 consecutive nodes |
| **Node drop** | Remove a customer if its insertion cost exceeds its profit |
| **Node add** | Insert an unvisited customer at its cheapest position if profitable |

First-improvement strategy: accept the first improving move found, restart the pass.

### Reduced graph

For the LP-guided callback, construction and local search operate on a **reduced graph** — a subset of the original problem's edges and nodes, selected based on LP relaxation values. The `edge_active` mask is threaded through `find_edge`, `edge_cost`, `greedy_insert`, and `local_search`; inactive edges return cost=infinity and are skipped.

Three reduction strategies are used concurrently:

**Strategy A: LP-value threshold**
- Include edges with `x_e > 0.1` + edges between pairs of nodes with `y_i > 0.5`
- Simple, focuses on the LP support subgraph

**Strategy B: RINS-style**
- Include edges from incumbent (`x_e > 0.5`) OR fractional LP (`x_e > 0.1`)
- Falls back to Strategy A when no incumbent exists
- Searches in the union of incumbent and LP relaxation neighborhoods

**Strategy C: Neighborhood expansion**
- Seed nodes: endpoints of edges with `x_e > 0.3`
- Include seed edges + all edges between pairs of active nodes
- Focuses on the connected neighborhood around the LP support

All strategies: source/target always active.

### Parallel restarts (TBB)

The initial heuristic supports two modes, controlled by the `deterministic` solver option (default: `true`):

**Deterministic mode** (default):
- Pre-builds all orderings upfront: 3 deterministic + random shuffles with seeds `0, 1, 2, ...`
- Runs all restarts via `tbb::parallel_for` over the fixed set
- Best solution selected by iterating results in index order (deterministic tiebreaking)
- Restart count: `clamp(num_nodes, 20, 200)`
- Identical results across runs on any hardware

**Opportunistic mode** (`--deterministic false`):
- Multiple TBB workers race with `std::random_device` seeds
- Workers run until a time budget expires: `min(500ms, num_nodes * 10ms)`
- May find better bounds on fast multi-core hardware
- Results vary between runs

## Standalone Usage

```cpp
#include "heuristic/primal_heuristic.h"

// Deterministic (default): fixed restart count, reproducible
auto warm_start = heuristic::build_initial_solution(problem_, num_restarts);

// Opportunistic: time-based, non-deterministic
auto warm_start = heuristic::build_initial_solution(problem_, num_restarts, budget_ms);

// warm_start.objective = travel_cost - collected_profit
// warm_start.col_values = solution vector (edges + nodes)
```

Python:
```python
from rcspp_bac import build_warm_start
warm = build_warm_start(prob, time_budget_ms=500.0)
```

### LP-guided callback (during solve)

Registered via `HiGHSBridge::install_heuristic_callback()`:

1. Separator callback caches LP relaxation values (`x_vals`, `y_vals`) under a mutex
2. `kCallbackMipUserSolution` callback reads cached LP values, builds reduced graphs, runs heuristic
3. If improved solution found: `data_in->user_has_solution = true; data_in->user_solution = result`
4. HiGHS validates via `solutionFeasible()` + our SEC feasibility check in `addIncumbent()`

The callback skips the initial `kExternalMipSolutionQueryOriginAfterSetup` query (pre-solve heuristic already ran). Rate-limited to once per 200 B&B nodes to avoid overhead on hard instances.

The solution vector has size `num_edges + num_nodes`:
- `col_values[e] = 1` for edges in the route
- `col_values[m + v] = 1` for visited nodes (including source/target)

## Testing

6 C++ tests (`[heuristic]` tag) and 5 Python tests verify:
- Tour produces a valid closed loop with correct degree (degree 2 at all visited nodes)
- Path produces a valid open path (degree 1 at source/target, degree 2 at intermediates)
- Total demand of visited nodes respects vehicle capacity
- Source/target y-variables are set to 1
- Edge count equals visited-node count minus 1 (path) or equals visited-node count (tour)
- Objective is finite

```bash
./build/rcspp_algo_tests [heuristic]
```

## HiGHS Integration

In `Model::solve()` (model.cpp), the warm start is fed to HiGHS:

```cpp
HighsSolution start;
start.value_valid = true;
start.col_value = std::move(warm.col_values);
highs.setSolution(start);
```

### CLI options

- `--heuristic_callback true/false` (default: true) — enable/disable LP-guided callback
- `--heuristic_budget_ms N` (default: 20) — time budget per callback invocation in ms
- `--heuristic_strategy N` (default: 0) — 0=all strategies, 1=LP-threshold, 2=RINS, 3=neighborhood

## Performance

Tested on SPPRCLIB instances (verified against known optimal values from Jepsen et al. 2014):

| Instance | Heuristic obj | Optimal | Gap | Heuristic time |
|---|---|---|---|---|
| B-n45-k6-54 | -70,029 | -74,278 | 5.7% | 0.45s |
| B-n50-k8-40 | -7,412 | -12,832 | 42% | 0.50s |
| B-n52-k7-15 | -63,789 | -74,998 | 15% | 0.50s |
| A-n63-k10-44 | -26,280 | -32,561 | 19% | 0.50s |
| A-n69-k9-42 | -35,613 | -43,290 | 18% | 0.50s |

Solver speedup on B-n45-k6-54: **13.6s -> 8.7s** (261 vs 548 nodes).

## Possible improvements

- **Synchronize + diversify**: After a batch of restarts, collect the k best tours,
  measure diversity (e.g. symmetric difference of visited nodes), recombine diverse
  parents to seed the next batch. Essentially a memetic algorithm.
- **Regret-based insertion**: Instead of cheapest-insert, use regret-2 (insert the
  customer whose best position advantage over second-best is largest). Standard in VRP.
- **Perturbation**: After local search converges, perturb (e.g. random segment removal +
  re-insertion) and re-optimize. Iterated local search.
- **Edge cost cache**: Precompute n*n distance matrix to avoid repeated `find_edge` scans.
  Currently ~0.1ms for find_edge on n=69, but would help for n>200.
