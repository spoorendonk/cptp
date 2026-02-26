# Algorithms and Techniques

## Problem

The solver targets a variant of the **Resource Constrained Shortest Path Problem (RCSPP)** with vertex profits and optional vertices:

- A vehicle starts at a **source** and ends at a **target** node
- Each customer has a **profit** (reward for visiting), a **demand** (capacity consumed), and a **cost** (travel distance)
- The vehicle has limited **capacity** Q
- The objective is to minimize `travel_cost - collected_profit`, i.e., find a route that balances short routes with high-value customers

Not all customers need to be visited -- the solver decides which subset maximizes net benefit.

The tour special case (source == target) corresponds to the **Capacitated Profitable Tour Problem (CPTP)** of Jepsen et al. (2014).

The solver supports two modes:
- **Tour** (source == target): closed loop from a depot, as in the original CPTP formulation
- **s-t path** (source != target): open path from source to target

## MIP Formulation

Undirected edge formulation following Jepsen et al. (2014):

- $x_e \in \{0,1\}$: whether edge $e$ is in the tour
- $y_i \in \{0,1\}$: whether customer $i$ is visited

**Objective**:

$$\min \sum_e c_e \, x_e - \sum_i p_i \, y_i$$

**Constraints**:
1. Degree:
   - **Tour**: $\sum_{e \in \delta(i)} x_e = 2\,y_i$ for each node $i$
   - **Path**: $\sum_{e \in \delta(i)} x_e = y_i$ for source/target, $= 2\,y_i$ for intermediates
2. Capacity: $\sum_i d_i \, y_i \le Q$
3. Fixed nodes: $y_s = 1$, $y_t = 1$
4. Subtour elimination (separated dynamically)
5. Capacity inequalities (separated dynamically)

## Solver Pipeline

```
Load instance
    |
    v
Preprocessing (parallel via TBB)
    |-- Forward labeling (bounds from source)
    |-- Backward labeling (bounds from target; same as forward for tour)
    '-- Initial heuristic (greedy + local search)
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
    |-- Feasibility check (SEC on incumbents, with sub-MIP column mapping)
    |-- Sub-MIP root SEC separation (column mapping via undoPrimal)
    |-- Domain propagator --> domain-propagator.md
    '-- LP-guided heuristic callback --> warm-start-heuristic.md
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
RGLM, and Comb. All cuts use formulations valid for RCSPP with optional vertices,
originally developed for the CPTP (Jepsen et al. 2014). A shared Gomory-Hu tree
avoids redundant max-flow computation.

See [separation.md](separation.md) for full details.

## Domain Propagator

A custom propagator injected into HiGHS via `HighsUserPropagator`. Uses labeling
bounds to fix edge variables to zero when the upper bound tightens or edges are
fixed by branching.

See [domain-propagator.md](domain-propagator.md) for full details.

## Preprocessing

Demand-reachability filtering and labeling-based edge elimination. See [preprocessing.md](preprocessing.md) for full details.

## Primal Heuristic

Parallel construction + local search heuristic using Intel TBB. See [warm-start-heuristic.md](warm-start-heuristic.md) for full details.

**Initial solution** (`build_initial_solution`): Greedy cheapest-insertion with three deterministic orderings (profit/demand ratio, absolute profit, distance from source) plus random restarts on the complete graph. Time budget: $\min(500\text{ms},\; n \times 10\text{ms})$.

**LP-guided callback** (`lp_guided_heuristic`): During MIP solve, builds reduced graphs from LP relaxation values using three strategies (LP-value threshold, RINS-style agreement, neighborhood expansion) and runs construction + local search on each. Registered via `kCallbackMipUserSolution`. Time budget: 20ms per invocation.

**Local search**: 2-opt, or-opt (1/2/3 chains), node drop, node add with first-improvement.

**Parallelism**: By default (deterministic mode), all restarts are pre-built with deterministic seeds and executed via `tbb::parallel_for`, giving identical results across runs. With `--deterministic false`, the solver uses time-based workers with random seeds for potentially better bounds on fast hardware.

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). A branch-and-cut algorithm for the capacitated profitable tour problem. *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Gusfield, D. (1990). Very simple methods for all pairs network flow analysis. *SIAM Journal on Computing*, 19(1), 143-155.
- SPPRCLIB: Benchmark instance library for shortest path problems with resource constraints. https://or.rwth-aachen.de/spprclib
- HiGHS: High-performance open-source linear and mixed-integer programming solver. https://highs.dev
