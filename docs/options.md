# Solver options

Every option is a `--key value` pair on `cptp-solve`, the same pair on
`python -m cptp`, and a `("key", "value")` tuple in Python's `Model.solve`
(where the source and target are set on the problem instead). The
block below is `cptp-solve --help`, word for word: a test fails when the two
differ, and `cmake --build build --target cptp_options_doc_update` rewrites it.

<!-- cptp-solve --help: begin -->
```text
cptp-solve — Capacitated Profitable Tour Problem solver

Usage: cptp-solve <instance> [--source <node>] [--target <node>] [--<option> <value> ...]

An instance has m >= 1 capacity rows (row 0 is the first), and a solution
satisfies every row; the options below say which rows each component uses.

Boolean options (shown as true/false) also take 1/0 or on/off, in any case;
any other value is an error.

General:
  --source <node>                          Source node (overrides file value)
  --target <node>                          Target node (overrides file value)
                                           If source != target, solves s-t path
  --cutoff <value>                         Known UB for prove-only mode
  --time_limit <sec>                       Wall-clock time limit (adjusted for preprocessing)
  --threads N                              HiGHS thread count
  --random_seed N                          Random seed for HiGHS and heuristics (default 0)
  --output_flag true/false                 HiGHS console output (default true)

Bounds & Preprocessing:
  --labeling_elim_bounds true/false        Capacity-aware labeling bounds (default true)
                                           When false, disables edge elim & propagation
  --labeling_max_queue_pops N              Max queue pops per labeling run
                                           (default 1000000, 0=unlimited)
  --labeling_rows first|all|surrogate|vector
                                           Capacity rows the labeling bounds use (default first:
                                           row 0; all: row 0 and every bounded row, combined;
                                           surrogate: row 0 and sum_k d_k/Q_k; vector: loads on every row)
  --labeling_join min|load                 Edge elimination joins the labeling halves at
                                           their cheapest labels (min) or by row-0 load (default load)
  --labeling_2cycle true/false             2-cycle elimination in the preprocessing single-row
                                           labelings (default false; not the vector or in-tree RC ones)
  --all_pairs_bounds true/false            All-pairs bounds for propagation + SPI
                                           (default false)
  --edge_elimination true/false            Edge elimination (default true, needs bounds)
  --edge_elimination_nodes true/false      Node fixing from eliminated edges
                                           (default true, needs bounds)

Propagation:
  --bounds_propagation true/false          Bounds-based domain propagation
                                           (default true, needs bounds)
  --rc_fixing <strategy>                   Reduced-cost fixing strategy
                                           off|root_only|on_ub_improvement|periodic|adaptive
  --rc_fixing_interval N                   Interval for periodic strategy (default 100)
  --rc_fixing_to_one true/false            Fix-to-1 for node variables
                                           (default false, expensive)

Heuristics:
  Warm-start (construction + local search):
  --heu_ws true/false                      Enable warm-start heuristic (default true)
  --heu_ws_ls_max_iter N                   Local search iteration limit (default 1000)

  LP-guided (B&C callback):
  --heu_lpg true/false                     Enable LP-guided heuristic (default true)
  --heu_lpg_strategy N                     0=all, 1=LP-threshold, 2=RINS, 3=neighborhood
  --heu_lpg_budget_ms N                    Time budget per callback in ms (default 20)
  --heu_lpg_node_interval N                Callback cadence in B&B nodes (default 200)
  --heu_lpg_deterministic_restarts N       Fixed restarts per callback
  --heu_lpg_edge_threshold X               LP-threshold edge cutoff (default 0.1)
  --heu_lpg_node_threshold X               LP-threshold node cutoff (default 0.5)
  --heu_lpg_lp_threshold X                 RINS LP edge cutoff (default 0.1)
  --heu_lpg_seed_threshold X               Neighborhood seed edge cutoff (default 0.3)

  HiGHS sub-MIP:
  --heu_highs_submip_sec true/false        SEC in HiGHS sub-MIPs (default true)

Separation (global):
  --separation_interval N                  Separate every N-th callback (default 1)
  --separation_tol X                       Fractional separation tolerance (default 0.1)
  --max_cuts_per_round N                   Global cut cap: pool all cuts, sort by
                                           violation, add top N (default 10, 0=unlimited)
  --root_sepa_rounds highs|unlimited|N     Round cap of HiGHS's root cut loop
                                           (default highs: 2*sqrt(log2 tree size); first root pass only)
                                           The loop also ends on a round without cuts, which
                                           separation_interval > 1 can cause
  --root_sepa_stall highs|unlimited|N      Stalled rounds before the root cut loop
                                           stops (default highs: 3)

Separation (per cut family):
  --enable_<family> true/false             Enable/disable a cut family

  --enable_sec true/false                  SEC (default true)
  --enable_rci true/false                  RCI (default true)
  --enable_multistar true/false            Multistar (default true)
  --enable_comb true/false                 Comb (default false)
  --enable_rglm true/false                 RGLM (default false)
  --capacity_cut_rows all|most_violated|first
                                           Capacity rows RCI, multistar and RGLM cuts come from (default all:
                                           every row with 0 < Q < 1e17; most_violated: per set, the
                                           row with the most violated cut; first: row 0)
  --capacity_cut_surrogate true/false      RCI and RGLM also cut on the surrogate row
                                           sum_k d^k/Q_k of the bounded rows (default false)
  --multistar_largest_side true/false      Multistar also cuts on each target's largest
                                           depot-target min-cut side (default false)
  --enable_xcover true/false               Cross-row lifted cover w-GLMs, covers lifted
                                           against every bounded row (default false)

  SPI (requires --all_pairs_bounds true):
  --enable_spi true/false                  SPI (default false)
  --spi_max_dp_size N                      Largest node set for exact Held-Karp
                                           lower bound; O(2^N·N²) (default 15)
  --spi_max_seeds N                        Infeasible pairs to grow into larger
                                           cuts, most-violated first (default 50)
  --spi_max_grow_candidates N              High-y nodes tried when extending a
                                           seed pair (default 40)

Branching:
  --branch_hyper <mode>                    Hyperplane branching mode
                                           off|pairs|clusters|demand|cardinality|all
                                           One integer column z = a.y per
                                           hyperplane, branched on by HiGHS

HiGHS:
  Unrecognized --key value pairs are forwarded to HiGHS.
  Presolve is always disabled (cptp requires stable column mapping).
  --highs_help                             Print all HiGHS options and exit
  --version                                Print cptp-solve and HiGHS versions, then exit
```
<!-- cptp-solve --help: end -->

## Defaults the help does not show

- `rc_fixing`: `adaptive`.
- `heu_lpg_strategy`: `0` (all three neighbourhoods).
- `heu_lpg_deterministic_restarts`: `32`.
- `time_limit`: no limit.
- `threads`: half the machine's hardware threads, rounded up, for HiGHS and
  for cptp's own parallel work alike.

## Options that go to HiGHS

A key cptp does not know goes to HiGHS, and one HiGHS rejects is logged as a
warning and ignored. Two that matter here:

- `parallel`: HiGHS's parallel tree search, on (HiGHS's `choose`) whenever a
  solve has more than one thread. With it, node counts, statistics and, among
  tied optima, the tour can differ between runs with the same seed;
  `parallel off` (or `threads 1`) gives a reproducible solve, with cptp's
  separators and labeling parallel inside one tree instead.
- `mip_rel_gap`, `mip_abs_gap`: stop at a gap instead of at a proof of
  optimality.

`presolve` cannot be turned on: cptp needs HiGHS's columns to stay as it built
them.

## Python

`python -m cptp` sets its own defaults: `time_limit` 600 s, `threads` 1, and
`output_flag` false (`--verbose` turns it on). `cptp.solve(...)` takes
`time_limit`, `num_threads` and `verbose` as arguments; any other option goes
through `Model.solve`.
