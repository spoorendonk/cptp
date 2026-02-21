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

## Separation (Cut Generation)

Three families of cuts are separated dynamically via HiGHS MIP callbacks. All separators share a **Gomory-Hu tree** built once per separation round, avoiding redundant max-flow computations.

### Subtour Elimination Constraints (SEC)

For each customer t, check if the min-cut between source and t in the fractional support graph violates the SEC. The Gomory-Hu tree (rooted at source) provides all pairwise min-cuts efficiently.

- **Tour**: `x(delta(S)) >= 2 * y_t` for any S not containing the depot
- **Path**: `x(delta(S)) >= 2 * y_t` if S does not contain the path target (path must enter and leave S); `x(delta(S)) >= y_t` if S contains the path target (path enters S and terminates)

Uses Dinitz max-flow to find violated cuts. When the sparser inside form has fewer nonzeros, a degree-substituted variant is used instead of the cut form.

### Rounded Capacity Inequalities (RCI)

For a set S of customers:

```
x(delta(S)) >= 2 * ceil(d(S) / Q)
```

where `d(S)` is the total demand of S. Strengthens the LP relaxation by accounting for the number of vehicle trips needed to serve S.

### Multistar Inequalities

Strengthens capacity cuts using knapsack-like reasoning on individual edge contributions. Uses demand-weighted coefficients to tighten the bound.

### Cut Management

- Cuts are sorted by violation magnitude; only the top-k (default 3) per separator per round are added (Jepsen et al. recommend keeping only the most-violated cuts)
- Default fractional separation tolerance is 0.1 (Jepsen et al. (2008) found 0.4 optimal with CPLEX's built-in cuts; HiGHS needs a lower threshold since it has fewer built-in cutting planes)

## Gomory-Hu Tree

Implements **Gusfield's simplified algorithm** to compute all pairwise min-cuts with only n-1 max-flow calls (instead of O(n^2)). The tree is built once per separation round and shared across all three separators.

Each max-flow call uses the **Dinitz algorithm** (O(V^2 E)) on the fractional support graph.

**References**:
- Gusfield, D. (1990). Very simple methods for all pairs network flow analysis. *SIAM Journal on Computing*, 19(1), 143-155.
- Dinitz, Y. (1970). Algorithm for solution of a problem of maximum flow in a network with power estimation. *Soviet Mathematics Doklady*, 11, 1277-1280.

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
