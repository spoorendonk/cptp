# rcspp-solve CLI

This document tracks RCSPP-specific options intercepted by `Model::solve()`
and CLI overrides handled in `src/cli/main.cpp`.

All options use `--key value`. Boolean values accept `true|false` (or `1|0`).

## Scope and forwarding

- RCSPP options listed here are intercepted by `rcspp-bac`.
- Any other `--key value` pair is forwarded to HiGHS.
- `presolve` is intercepted and always forced to `off`.

## Invocation and runtime controls

| Option | Default | Values | Effect |
|---|---:|---|---|
| `source` | from instance | node id | CLI-only override for source node (`--source`). |
| `target` | from instance | node id | CLI-only override for target node (`--target`). |
| `cutoff` | unset | finite double | Sets known UB/prove-only bound (`objective_bound` in HiGHS). |
| `output_flag` | `true` | `true/false` | Enables/disables solver log forwarding to rcspp logger. |
| `presolve` | forced off | any | Any user request is ignored; solver runs with `presolve=off`. |

Notes:
- `--source/--target` are parsed by `main.cpp`, not by `Model::solve()` options.
- The solver internally sets HiGHS `output_flag=false` and forwards logs through
  the logger only when `output_flag=true`.

## 1) Cut-separation policy

| Option | Default | Values | Effect |
|---|---:|---|---|
| `enable_sec` | `true` | `true/false` | Enable SEC cut generation in main MIP and sub-MIP separation paths. |
| `enable_rci` | `true` | `true/false` | Enable RCI separator. |
| `enable_multistar` | `true` | `true/false` | Enable Multistar separator. |
| `enable_comb` | `true` | `true/false` | Enable Comb separator. |
| `enable_spi` | `true` | `true/false` | Enable SPI separator when `all_pairs_propagation=true`. |
| `max_cuts_sec` | unset | integer (`>=0`) | Per-round cap for SEC cuts. |
| `max_cuts_rci` | unset | integer (`>=0`) | Per-round cap for RCI cuts. |
| `max_cuts_multistar` | unset | integer (`>=0`) | Per-round cap for Multistar cuts. |
| `max_cuts_comb` | unset | integer (`>=0`) | Per-round cap for Comb cuts. |
| `max_cuts_rglm` | unset | integer (`>=0`) | Per-round cap for RGLM cuts. |
| `max_cuts_spi` | unset | integer (`>=0`) | Per-round cap for SPI cuts. |
| `min_violation_sec` | unset | double (`>=0`) | Minimum SEC violation required before adding cut. |
| `min_violation_rci` | unset | double (`>=0`) | Minimum RCI violation required before adding cut. |
| `min_violation_multistar` | unset | double (`>=0`) | Minimum Multistar violation required before adding cut. |
| `min_violation_comb` | unset | double (`>=0`) | Minimum Comb violation required before adding cut. |
| `min_violation_rglm` | unset | double (`>=0`) | Minimum RGLM violation required before adding cut. |
| `min_violation_spi` | unset | double (`>=0`) | Minimum SPI violation required before adding cut. |
| `enable_rglm` | `false` | `true/false` | Enable RGLM separator. |
| `separation_interval` | `1` | integer (`>=1` effective) | Run separation every Nth separator callback. |
| `max_cuts_per_separator` | `3` | integer (`>=0`) | Global per-separator cut cap (`0` means unlimited). |
| `separation_tol` | `0.1` | nonnegative double | Fractional violation threshold. |
| `submip_separation` | `true` | `true/false` | Enable SEC separation in sub-MIP root. |

Notes:
- Per-family caps fall back to `max_cuts_per_separator`.
- Per-family min-violation thresholds fall back to `0.0`.
- Even with `enable_sec=false`, SEC feasibility checks still protect incumbent acceptance.

See also: [Cut Separation](separation.md).

## 2) Heuristics (warm start + in-search)

| Option | Default | Values | Effect |
|---|---:|---|---|
| `disable_heuristics` | `false` | `true/false` | Disable warm-start injection and set `mip_heuristic_effort=0`. |
| `deterministic` | `true` | `true/false` | Deterministic vs opportunistic preprocessing/heuristic behavior. |
| `heuristic_callback` | `true` | `true/false` | Enable LP-guided callback heuristic. |
| `heuristic_budget_ms` | `20.0` | positive double | Callback time budget (ms). |
| `heuristic_strategy` | `0` | `0/1/2/3` | `0=all, 1=LP-threshold, 2=RINS, 3=neighborhood`. |
| `heuristic_node_interval` | `200` | integer | Callback cadence in B&B nodes (clamped to `>=1`). |
| `heuristic_async_injection` | `true` | `true/false` | Allow async incumbent injection from background updates. |

See also: [Primal Heuristic](primal-heuristic.md).

## 3) Propagation, ng bounds, and preprocessing

| Option | Default | Values | Effect |
|---|---:|---|---|
| `edge_elimination` | `true` | `true/false` | Enable UB-based edge elimination preprocessing. |
| `edge_elimination_nodes` | `true` | `true/false` | Also fix node vars when all incident edges are eliminated. |
| `all_pairs_propagation` | `false` | `true/false` | Enable all-pairs Trigger B propagation (and SPI availability). |
| `dssr_async` | `true` | `true/false` | Background DSSR bound tightening during B&B. |
| `ng_initial_size` | `4` | integer | Initial ng size (clamped to `>=1`). |
| `ng_max_size` | `12` | integer | Max ng size (clamped to `>= ng_initial_size`). |
| `ng_dssr_iters` | `6` | integer | DSSR iterations (clamped to `>=1`). |
| `ng_label_budget` | `50` | integer | Max active labels per node (clamped to `>=1`). |
| `ng_simd` | `true` | `true/false` | Enable AVX2 prefilter path in labeling dominance scans. |

Notes:
- Async DSSR worker is only launched when `dssr_async=true` and
  `all_pairs_propagation=false`.

See also: [Preprocessing](preprocessing.md), [Domain Propagator](domain-propagator.md), [ng/DSSR Bounds](ng-dssr-labelling.md).

## 4) Reduced-cost fixing

| Option | Default | Values | Effect |
|---|---:|---|---|
| `rc_fixing` | `adaptive` | `off/root_only/on_ub_improvement/periodic/adaptive` | Reduced-cost fixing trigger policy. |
| `rc_fixing_interval` | `100` | integer | Period for `periodic` mode (use `>=1`). |
| `rc_fixing_to_one` | `false` | `true/false` | Enable expensive node fix-to-1 pass. |

See also: [Domain Propagator](domain-propagator.md).

## 5) Hyperplane branching

| Option | Default | Values | Effect |
|---|---:|---|---|
| `branch_hyper` | `off` | `off/pairs/clusters/demand/cardinality/all` | Select hyperplane candidate family. |
| `branch_hyper_sb_max_depth` | `0` | integer | Depth cap for hyperplane strong-branch probing (clamped to `>=0`). |
| `branch_hyper_sb_iter_limit` | `100` | integer | LP iteration cap per strong-branch trial (clamped to `>=1`). |
| `branch_hyper_sb_min_reliable` | `4` | integer | Pseudocost reliability threshold (clamped to `>=1`). |
| `branch_hyper_sb_max_candidates` | `3` | integer | Max strong-branch candidate count (clamped to `>=1`). |

See also: [Hyperplane Branching](hyperplane-branching.md).

## 6) Concurrent solve admission

| Option | Default | Values | Effect |
|---|---:|---|---|
| `max_concurrent_solves` | `0` | integer | In-process concurrent solve cap (`0` means no per-call cap; clamped to `>=0`). |

See also: [Solve Concurrency](solve-concurrency.md).

## Example bundles

### A. Conservative/stable debug profile

```bash
--deterministic true \
--max_concurrent_solves 1 \
--enable_sec true --enable_rci false --enable_multistar false --enable_comb false --enable_spi false \
--max_cuts_sec 1 --min_violation_sec 1e-4
```

### B. Strong cut-separation profile

```bash
--enable_sec true --enable_rci true --enable_multistar true --enable_comb true --enable_spi true \
--max_cuts_sec 5 --max_cuts_rci 5 --max_cuts_multistar 5 --max_cuts_comb 5 --max_cuts_spi 5 \
--min_violation_sec 1e-4 --min_violation_rci 1e-4 --min_violation_multistar 1e-4 --min_violation_comb 1e-4 --min_violation_spi 1e-4
```

### C. Propagation-heavy profile

```bash
--all_pairs_propagation true \
--edge_elimination true --edge_elimination_nodes true \
--ng_initial_size 4 --ng_max_size 16 --ng_dssr_iters 8 --ng_label_budget 100 --dssr_async true
```

### D. Reduced-cost fixing profile

```bash
--rc_fixing periodic --rc_fixing_interval 50 --rc_fixing_to_one false
```

### E. Heuristic-light prove-only profile

```bash
--cutoff <known_ub> \
--disable_heuristics true \
--heuristic_callback false
```
