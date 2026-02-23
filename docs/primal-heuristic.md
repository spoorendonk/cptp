# Primal Heuristic

File: `src/heuristic/primal_heuristic.h` (header-only)

## Purpose

Two heuristic components:
1. **Initial solution** (`build_initial_solution`): runs before MIP solve on the complete graph, provides an initial feasible solution to HiGHS via `setSolution()`
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
- Include edges with `x_e > 0.1` and all edges incident to nodes with `y_i > 0.5`
- Simple, preserves edges the LP uses + neighborhood around visited nodes

**Strategy B: RINS-style agreement**
- Include edges where LP and incumbent agree: both `> 0.5` or both `< 0.5`
- Falls back to Strategy A when no incumbent exists
- Searches near where LP and best solution agree

**Strategy C: Neighborhood expansion**
- Seed: edges with `x_e > 0.3`
- Expand 1-hop: add all edges incident to any node touched by seed edges
- Gives a connected neighborhood around the LP support

All strategies: source/target always active.

### Parallel restarts (TBB)

Multiple workers run construction + local search concurrently:
- Each worker picks the next strategy (round-robin) and ordering
- First round per strategy: LP-weighted ordering (sort by `y_lp` descending, break ties by profit/demand ratio)
- Subsequent rounds: random shuffles within each strategy
- All workers share a best-known solution (mutex-protected)
- Workers run until a time budget expires

Initial heuristic budget: `min(500ms, num_nodes * 10ms)`.
LP-guided callback budget: 20ms per invocation (configurable).

## Integration

### Initial solution (pre-solve)

In `Model::solve()` (model.cpp):

```cpp
auto warm_start = heuristic::build_initial_solution(problem_, budget_ms);
HighsSolution start;
start.value_valid = true;
start.col_value = std::move(warm_start.col_values);
highs.setSolution(start);
```

### LP-guided callback (during solve)

Registered via `HiGHSBridge::install_heuristic_callback()`:

1. Separator callback caches LP relaxation values (`x_vals`, `y_vals`) under a mutex
2. `kCallbackMipUserSolution` callback reads cached LP values, builds reduced graphs, runs heuristic
3. If improved solution found: `data_in->user_has_solution = true; data_in->user_solution = result`
4. HiGHS validates via `solutionFeasible()` + our SEC feasibility check in `addIncumbent()`

The callback skips the initial `kExternalMipSolutionQueryOriginAfterSetup` query (pre-solve heuristic already ran).

The solution vector has size `num_edges + num_nodes`:
- `sol[e] = 1` for edges in the route
- `sol[m + v] = 1` for visited nodes (including source/target)

### CLI options

- `--heuristic_callback true/false` (default: true) — enable/disable LP-guided callback
- `--heuristic_budget_ms N` (default: 20) — time budget per callback invocation in ms

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
