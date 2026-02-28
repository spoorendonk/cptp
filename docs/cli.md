# rcspp-solve CLI (Future Reference)

This document tracks the RCSPP-specific tuning options currently intercepted by `Model::solve` and intended for CLI exposure.

All options use `--key value` form. Boolean values should be passed as `true|false` (or `1|0`).

## Scope and forwarding

- RCSPP options listed here are intercepted and handled by `rcspp-bac`.
- Any other `--key value` pair is forwarded to HiGHS.
- `presolve` is a special case: it is intercepted and forced to `off` (details below).

## 1) Cut-separation policy

Use when you want to enable/disable cut families or tune cut aggressiveness per family.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `enable_sec` | `true` | `true/false` | Enable SEC cut generation in main MIP and sub-MIP separation paths. |
| `enable_rci` | `true` | `true/false` | Enable RCI separator (tour mode only). |
| `enable_multistar` | `true` | `true/false` | Enable Multistar separator (tour mode only). |
| `enable_comb` | `true` | `true/false` | Enable Comb separator (tour mode only). |
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

Notes:
- If a per-family cap is unset, it falls back to `max_cuts_per_separator`.
- If a per-family min violation is unset, it falls back to `0.0` for that family.
- Even with `enable_sec=false`, incumbent feasibility still checks SEC connectivity before accepting solutions.
- In `s-t` path mode, `RCI/Multistar/Comb/RGLM` are not added (tour-focused families).
- `SPI` is only active when `all_pairs_propagation=true`.

Related existing knobs:
- `enable_rglm` (default `false`)
- `separation_interval` (default `1`)
- `max_cuts_per_separator` (default `3`)
- `separation_tol` (default `0.1`)
- `submip_separation` (default `true`)

## 2) Heuristics (warm start + in-search)

Use when balancing startup quality, determinism, and heuristic overhead inside B&B.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `disable_heuristics` | `false` | `true/false` | Disables initial warm-start injection and disables HiGHS built-in MIP heuristics (`mip_heuristic_effort=0`). |
| `parallel_mode` | `deterministic` | `deterministic/opportunistic` | Controls solver-side parallel behavior policy. `deterministic`: fixed-restart deterministic heuristics/setup. `opportunistic`: randomized/time-budgeted behavior. |
| `deterministic_work_units` | `0` | integer (`>=0`) | Global non-HiGHS work budget for warm-start, callback heuristics, and async DSSR stages. `0` = uncapped count-only. |
| `heuristic_deterministic_restarts` | `32` | integer (`>=1`) | Fixed restart count per LP-guided callback invocation in deterministic mode. |
| `heuristic_node_interval` | `200` | integer (`>=1`) | Run LP-guided heuristic callback every N B&B nodes. |
| `heuristic_async_injection` | `true` | `true/false` | Allow async incumbent injection from background improvements in callback. |

Related existing knobs:
- `heuristic_callback` (default `true`)
- `heuristic_budget_ms` (default `20.0`)
- `heuristic_strategy` (`0=all, 1=LP-threshold, 2=RINS, 3=neighborhood`)

Notes:
- `deterministic_work_units` counts and optionally caps non-HiGHS activity only; HiGHS internal work is still controlled by HiGHS options (`threads`, limits, etc.).
- In deterministic mode, callback heuristics switch from wall-clock budgets to fixed restart counts.
- Legacy boolean mode flags `deterministic` and `opportunistic` are removed; use `parallel_mode`.

## 3) Propagation, NG bounds, and preprocessing

Use when tuning bound-tightening strength vs overhead.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `edge_elimination` | `true` | `true/false` | Enables UB-based edge elimination preprocessing before solve. |
| `edge_elimination_nodes` | `true` | `true/false` | Also fixes node vars when all incident edges are eliminated. |

Related existing knobs in the same area:
- `all_pairs_propagation` (default `false`)
- `dssr_background_updates` (default `false` when `parallel_mode=deterministic`, `true` otherwise)
- `dssr_background_policy` (default `fixed`): `fixed` runs all planned async stages (subject to caps/budget); `auto` stops background DSSR early when there is no observed search/checkpoint activity or no recent DSSR improvements.
- `dssr_background_max_epochs` (default `0`): hard cap on async DSSR stages (`0` = uncapped)
- `dssr_background_auto_min_epochs` (default `4`): minimum async DSSR stages before `auto` early-stop checks
- `dssr_background_auto_no_progress_limit` (default `6`): `auto` stops after this many non-improving stages
- `ng_initial_size` (default `4`)
- `ng_max_size` (default `12`)
- `ng_dssr_iters` (default `6`)
- `ng_label_budget` (default `50`)
- `ng_simd` (default `true`)

Note:
- Edge elimination requires a finite incumbent upper bound plus labeling bounds; if no UB is available yet, elimination is skipped.

## 4) Reduced-cost fixing

Use when tightening domains from LP reduced costs during propagation.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `rc_fixing` | `adaptive` | `off/root_only/on_ub_improvement/periodic/adaptive` | Selects reduced-cost fixing trigger policy. |
| `rc_fixing_interval` | `100` | integer (`>=1` recommended) | Period for `periodic` mode. |
| `rc_fixing_to_one` | `false` | `true/false` | Also tries fix-to-1 on node variables (more expensive). |

## 5) Hyperplane branching strong-branch tuning

Use when `branch_hyper != off` and you want to change strong-branch effort.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `branch_hyper_sb_max_depth` | `0` | integer (`>=0`) | Depth cap for hyperplane strong-branch probing. |
| `branch_hyper_sb_iter_limit` | `100` | integer (`>=1`) | LP iteration limit per hyperplane strong-branch trial. |
| `branch_hyper_sb_min_reliable` | `4` | integer (`>=1`) | Pseudocost reliability threshold before reducing strong branching. |
| `branch_hyper_sb_max_candidates` | `3` | integer (`>=1`) | Max hyperplane candidates tested with strong branching. |

Related existing knob:
- `branch_hyper`: `off/pairs/clusters/demand/cardinality/all`

## 6) Concurrent solve admission

Use when multiple threads/process components call `Model::solve` concurrently in the same process.

| Option | Default | Values | Effect |
|---|---:|---|---|
| `max_concurrent_solves` | `0` | integer (`>=0`) | Per-call cap for in-process solve admission. `0` means this caller adds no cap. |

Behavior:
- A stricter cap requested by any active or waiting solve call constrains later arrivals.
- Effective admission behaves like a shared strict-cap gate (minimum active/pending finite cap).
- This is independent from HiGHS `threads` (internal parallelism of one solve).

## 7) Presolve capture behavior

| Option | Effective behavior |
|---|---|
| `presolve` | Captured by RCSPP layer. Any request to enable presolve is ignored and a warning is logged. Solve always runs with HiGHS `presolve=off`. |

This is intentional for now because callback/cut plumbing assumes stable original row/column mappings.

## Example bundles

### A. Conservative/stable debug profile

```bash
--parallel_mode deterministic \
--deterministic_work_units 500 \
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
--parallel_mode opportunistic \
--all_pairs_propagation true \
--edge_elimination true --edge_elimination_nodes true \
--ng_initial_size 4 --ng_max_size 16 --ng_dssr_iters 8 --ng_label_budget 100 --dssr_background_updates true
```

### C2. Easy-instance guarded async DSSR profile

```bash
--parallel_mode opportunistic \
--dssr_background_updates true \
--dssr_background_policy auto \
--dssr_background_max_epochs 8
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
