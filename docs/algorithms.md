# Algorithms and Techniques

## Problem

The **Capacitated Profitable Trip Problem (CPTP)** is a variant of the Traveling Salesman Problem where:

- A vehicle starts at a **source** and ends at a **target** node
- Each customer has a **profit** (reward for visiting), a **demand** (capacity consumed), and a **cost** (travel distance)
- The vehicle has limited **capacity** Q
- The objective is to minimize `travel_cost - collected_profit`, i.e., find a route that balances short routes with high-value customers

Not all customers need to be visited -- the solver decides which subset maximizes net benefit.

The solver supports two modes:
- **Tour** (source == target): closed loop from a depot, as in the original CPTP formulation
- **s-t path** (source != target): open path from source to target

## MIP Formulation

Undirected edge formulation following Jepsen et al. (2014):

- **x_e in {0,1}**: whether edge e is in the tour
- **y_i in {0,1}**: whether customer i is visited

**Objective**: `min sum(cost_e * x_e) - sum(profit_i * y_i)`

**Constraints**:
1. Degree:
   - **Tour**: `sum(x_e : e incident to i) = 2 * y_i` for each node
   - **Path**: `sum(x_e : e incident to i) = y_i` for source/target, `= 2 * y_i` for intermediates
2. Capacity: `sum(demand_i * y_i) <= Q`
3. Fixed nodes: `y_source = 1`, `y_target = 1`
4. Subtour elimination (separated dynamically)
5. Capacity inequalities (separated dynamically)

## Solver Pipeline

```
Load instance
    |
    v
Preprocessing (parallel via TBB)
    |-- Forward labeling (ESPPRC bounds from source)
    |-- Backward labeling (bounds from target; same as forward for tour)
    '-- Warm-start heuristic (greedy + local search)
    |
    v
Build MIP formulation
    |-- Edge variables x_e, node variables y_i
    |-- Degree constraints, capacity constraint
    '-- Fix eliminated edges to zero (using labeling bounds + warm-start UB)
    |
    v
Install HiGHS callbacks
    |-- User separator   --> separation.md
    |-- Feasibility check (SEC on incumbents)
    '-- Domain propagator --> domain-propagator.md
    |
    v
HiGHS MIP solve
    |
    v
Extract result (tour/path, objective, bound, gap)
```

Orchestrated by `Model::solve()` in `src/model/model.cpp`.

## Cut Separation

Five families of cutting planes are separated dynamically: SEC, RCI, Multistar/GLM,
RGLM, and Comb. All cuts use CPTP-valid formulations (Jepsen et al. 2014) that
account for optional nodes. A shared Gomory-Hu tree avoids redundant max-flow
computation.

See [separation.md](separation.md) for full details.

## Domain Propagator

A custom propagator injected into HiGHS via `HighsUserPropagator`. Uses ESPPRC
labeling bounds to fix edge variables to zero when the upper bound tightens or
edges are fixed by branching.

See [domain-propagator.md](domain-propagator.md) for full details.

## Preprocessing

### Demand-Reachability

- **Tour**: Dijkstra from depot; eliminates nodes where `2 * min_demand_path(depot, v) - demand(v) > Q` (round-trip check)
- **Path**: Dijkstra from both source and target; eliminates nodes where `dist_source(v) + dist_target(v) - demand(v) > Q` (one-way check)

### Edge Elimination

ESPPRC-style forward labeling with:
- Net cost tracking (edge costs minus collected profits)
- Demand accumulation with capacity check
- 2-cycle elimination (no immediate return to predecessor)
- Label dominance to control enumeration

- **Tour**: labeling from depot; lb(u,v) = f[u] + c(u,v) + f[v] + profit(depot)
- **Path**: labeling from both source and target; lb(u,v) = min(f_s[u] + c(u,v) + f_t[v], f_s[v] + c(u,v) + f_t[u])

For each edge, if the lower bound exceeds the warm-start upper bound, the edge variable is fixed to zero.

## Warm-Start Heuristic

Parallel construction + local search heuristic using Intel TBB. See [warm-start-heuristic.md](warm-start-heuristic.md) for full details.

**Construction**: Greedy cheapest-insertion with three deterministic orderings (profit/demand ratio, absolute profit, distance from source) plus random restarts. For paths, the initial route is `[source, target]` (open); for tours, `[depot, depot]` (closed).

**Local search**: 2-opt, or-opt (1/2/3 chains), node drop, node add with first-improvement.

**Parallelism**: Multiple TBB workers run independent construction+search with a shared best solution. Time budget: `min(500ms, num_nodes * 10ms)`.

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). A branch-and-cut algorithm for the capacitated profitable tour problem. *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Gusfield, D. (1990). Very simple methods for all pairs network flow analysis. *SIAM Journal on Computing*, 19(1), 143-155.
- SPPRCLIB: Benchmark instance library for shortest path problems with resource constraints. https://or.rwth-aachen.de/spprclib
- HiGHS: High-performance open-source linear and mixed-integer programming solver. https://highs.dev
- PathWyse: Exact solver for resource-constrained shortest path problems (used for solution verification).
