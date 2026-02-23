# Preprocessing

Files: `src/preprocess/demand_reachability.h`, `src/preprocess/edge_elimination.h`

## Overview

Two preprocessing phases run in parallel (via TBB) before the MIP is built. Both
reduce problem size by eliminating nodes and edges that cannot appear in any
improving solution.

## Demand-Reachability

File: `src/preprocess/demand_reachability.h`

Eliminates nodes that are unreachable within the capacity budget.

- **Tour** (source == target): Dijkstra from depot; eliminates nodes where
  `2 * min_demand_path(depot, v) - demand(v) > Q` (round-trip check)
- **Path** (source != target): Dijkstra from both source and target; eliminates
  nodes where `dist_source(v) + dist_target(v) - demand(v) > Q` (one-way check)

## Edge Elimination

File: `src/preprocess/edge_elimination.h`

Forward labeling with:
- Net cost tracking (edge costs minus collected profits)
- Demand accumulation with capacity check
- 2-cycle elimination (no immediate return to predecessor)
- Label dominance to control enumeration

Lower bound computation:

- **Tour**: labeling from depot; `lb(u,v) = f[u] + c(u,v) + f[v] + profit(depot)`
- **Path**: labeling from both source and target;
  `lb(u,v) = min(f_s[u] + c(u,v) + f_t[v], f_s[v] + c(u,v) + f_t[u])`

For each edge, if the lower bound exceeds the warm-start upper bound, the edge
variable is fixed to zero before the MIP is built. The same labeling bounds are
reused by the [domain propagator](domain-propagator.md) during branch-and-bound.
