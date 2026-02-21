# Warm-Start Heuristic

File: `src/heuristic/warm_start.h` (header-only)

## Purpose

Provides an initial feasible solution to HiGHS via `setSolution()` before MIP solving.
A good primal bound prunes more of the B&B tree early.

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

### Parallel restarts (TBB)

Multiple workers run construction + local search concurrently:
- Each worker picks the next ordering (deterministic first, then random shuffles)
- All workers share a best-known solution (mutex-protected)
- Workers run until a time budget expires

Time budget: `min(500ms, num_nodes * 10ms)`
- 4 nodes: 40ms
- 45 nodes: 450ms
- 50+ nodes: 500ms

## Integration

In `Model::solve()` (model.cpp):

```cpp
auto warm_start = heuristic::build_warm_start(problem_, budget_ms);
// warm_start.objective = travel_cost - collected_profit
// warm_start.col_values = solution vector (edges + nodes)

HighsSolution start;
start.value_valid = true;
start.col_value = std::move(warm_start.col_values);
highs.setSolution(start);
```

The solution vector has size `num_edges + num_nodes`:
- `sol[e] = 1` for edges in the route
- `sol[m + v] = 1` for visited nodes (including source/target)

## Performance

Tested on SPPRCLIB instances (verified optimal against PathWyse):

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
