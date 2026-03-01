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
Stage 1 preprocessing (parallel via TBB)
    |-- Forward 2-cycle labeling from source  -> fwd bounds
    |-- Backward 2-cycle labeling from target -> bwd bounds (or fwd for tour)
    '-- Fast warm-start (quick UB, no bounds required)
    |
    v
First edge-elimination pass
    '-- Uses (Stage-1 bounds, fast UB)
    |
    v
Stage 2 preprocessing (parallel via TBB)
    |-- Second warm-start on reduced graph (adaptive, optional)
    |-- s-t ng/DSSR refinement (optional)
    '-- All-pairs 2-cycle bounds (optional)
    |
    v
Select final preprocessing artifacts
    |-- Final bounds snapshot (prefer refined ng, else all-pairs rows, else Stage-1)
    '-- Initial incumbent + UB
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
    '-- LP-guided heuristic callback --> primal-heuristic.md
    |
    v
Async ng/DSSR worker (optional, runs during HiGHS)
    |-- Publishes tighter (fwd,bwd,ng_size) snapshots
    '-- Publishes improved incumbents when found
    |
    v
HiGHS MIP solve
    |
    v
Extract result (tour/path, objective, bound, gap)
```

Orchestrated by `Model::solve()` in `src/model/model.cpp`. Set
`workflow_dump=true` to print startup/solve DAG wiring in logs.

## High-Level Pseudocode

```text
function solve(problem, options):
    configure_highs(options)

    # Stage 1: bounds + fast UB in parallel
    stage1_mode = choose(preproc_stage1_bounds in {two_cycle, ng1, auto})
    parallel:
        if stage1_mode == two_cycle:
            fwd = labeling_from(problem, source)
            bwd = labeling_from(problem, target)  # bwd=fwd in tour mode
        else:
            fwd, bwd, ng_path = ng_compute_bounds(problem, source, target,
                                                  ng_size=1, dssr_iters=1)
        ws_fast = build_initial_solution(problem, fast_budget)

    if source == target:
        bwd = fwd
    ub_fast = min(cutoff_if_any, ws_fast.objective)

    # First elimination estimate from (LB, UB)
    eliminated_stage1 = edge_elimination(problem, fwd, bwd, ub_fast, correction)
    run_second_ws = adaptive_decision(eliminated_stage1, n, ub_fast, options)

    # Stage 2: refinement and optional second warm-start
    parallel:
        if run_second_ws:
            ws_second = build_initial_solution(problem, full_budget, fwd, bwd,
                                              correction, ub_fast)
        if options.all_pairs_propagation:
            all_pairs = all_pairs_forward_labeling(problem)
        if ng_refinement_enabled(options):
            ng_bounds, ng_path = ng_compute_bounds(problem, source, target, options)

    (fwd_final, bwd_final, ng_size_used) = choose_bounds_snapshot(
        ng_bounds, all_pairs, fwd, bwd)
    warm_start = best_of(ws_fast, ws_second, ng_path_start_if_any)
    warm_start_ub = best_objective(cutoff_if_any, warm_start, ng_path)

    build_formulation_with_static_elimination(fwd_final, bwd_final, warm_start_ub)
    install_callbacks(separator, propagator, lp_heuristic)
    set_initial_solution(warm_start)

    if dssr_async and not all_pairs_propagation and ng_size_used < ng_max:
        launch_async_ng_worker(shared_bounds_store, async_incumbent_store)

    highs.run()
    stop_async_ng_worker()
    return extract_result()
```

## ParaMIP Worker DAG (Current `static_root` Direction)

`paramip_mode=static_root` executes a root-partitioned ParaMIP-style solve in-process.
The current flow is:

```text
stage1+stage2 preprocessing (once in parent)
      -> stage0_root_probes(seed_0..seed_k, mip_max_nodes=1, threads=1)
      -> pick_root_candidate(best|first)
      -> extract selected root LP (y-fractionality + basis + solution)
      -> create_branch_chunks(B) from most-fractional root y-vars
        -> worker_1 solve chunk_1 (threads=1)
        -> worker_2 solve chunk_2 (threads=1)
        -> ...
        -> worker_B solve chunk_B (threads=1)
          (each chunk reuses parent preprocessing and gets root basis/solution warm-start)
      -> aggregate(best incumbent + global bound/stats)
```

Current implementation status:
- `paramip_mode=plan`: builds and logs deterministic root chunk plans.
- `paramip_mode=static_root`: executes static root chunk solves and aggregates results.
- Stage-0 root probe races run with different random seeds and `mip_max_nodes=1`.
- Root probe winner policy: `paramip_root_pick=auto|best|first` (`auto` = best in deterministic, first in opportunistic).
- Chunk split variables are selected after Stage-0 from the selected root LP (`y` most-fractional first, fallback heuristic if needed).
- `parallel_mode=deterministic`: executes chunk solves in fixed chunk order.
- `parallel_mode=opportunistic`: executes chunk solves with a worker pool (`paramip_workers`).

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

Parallel construction + local search heuristic using Intel TBB. See [primal-heuristic.md](primal-heuristic.md) for full details.

**Initial solution** (`build_initial_solution`): Greedy cheapest-insertion with three deterministic orderings (profit/demand ratio, absolute profit, distance from source) plus random restarts. It is used in two startup roles: a fast first pass to get an early UB, and an optional second pass on a bound-filtered reduced graph after the first elimination step.

**LP-guided callback** (`lp_guided_heuristic`): During MIP solve, builds reduced graphs from LP relaxation values using three strategies (LP-value threshold, RINS-style agreement, neighborhood expansion) and runs construction + local search on each. Registered via `kCallbackMipUserSolution`. Time budget: 20ms per invocation.

**Local search**: 2-opt, or-opt (1/2/3 chains), node drop, node add with first-improvement.

**Parallelism**: By default (`--parallel_mode deterministic`), all restarts are pre-built with deterministic seeds and executed via `tbb::parallel_for`, giving identical results across runs. With `--parallel_mode opportunistic`, the solver uses time-based workers with random seeds for potentially better bounds on fast hardware.

## Hyperplane Branching

Optional hyperplane candidate branching over `y` variables is available via
`--branch_hyper` (`off/pairs/clusters/demand/cardinality/all`) with
strong-branch tuning parameters `branch_hyper_sb_*`.

See [hyperplane-branching.md](hyperplane-branching.md) for full details.

## Solve Concurrency

In-process solve admission can be controlled with
`--max_concurrent_solves`, independent from HiGHS thread count.

See [solve-concurrency.md](solve-concurrency.md) for full details.

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). A branch-and-cut algorithm for the capacitated profitable tour problem. *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Gusfield, D. (1990). Very simple methods for all pairs network flow analysis. *SIAM Journal on Computing*, 19(1), 143-155.
- SPPRCLIB: Benchmark instance library for shortest path problems with resource constraints. https://or.rwth-aachen.de/spprclib
- HiGHS: High-performance open-source linear and mixed-integer programming solver. https://highs.dev
