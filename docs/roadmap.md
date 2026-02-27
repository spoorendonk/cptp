# rcspp-bac — Implementation Roadmap

Branch-and-cut solver for the Resource Constrained Shortest Path Problem (RCSPP).

Current state: production-quality solver with SEC, RCI, Multistar, RGLM, Comb
separators, demand-reachability preprocessing, edge elimination via ng-labeling,
domain propagation, and parallel primal heuristic (warm-start). Solves 45/45 SPPRCLIB
instances and 26/31 Roberti Set 3 instances to optimality.

Latest update (2026-02-27): Step 0.5 and 0.6 completed — adaptive RC-fixing
policy study + defaults, and async incumbent/proof handoff from ng/DSSR into
HiGHS callbacks.

---

## 1. Current Architecture

```
src/core/        — Problem, static_graph, IO, Dinitz max-flow, Gomory-Hu tree
src/preprocess/  — Demand-reachability, ng/DSSR labeling, edge elimination
src/sep/         — SEC, RCI, Multistar, RGLM, Comb separators
src/model/       — HiGHS integration (Model, HiGHSBridge, propagator)
src/heuristic/   — Primal heuristic (warm-start) construction + local search
src/cli/         — CLI tool (rcspp-solve)
src/util/        — Logger, Timer
python/          — nanobind Python bindings
tests/           — Catch2 tests
docs/            — Algorithm documentation
```

Key design: all separators share a Gomory-Hu tree computed once per round
(Gusfield's algorithm, O(n) max-flows). Separators run in parallel via TBB
task_group. Domain propagator uses labeling bounds for variable fixing.

---

## 2. Milestones (History-Based Status)

1. [x] Core branch-and-cut solver foundation
2. [x] Demand-reachability preprocessing and amortized separation
3. [x] Edge elimination preprocessing via labeling
4. [x] Full separator stack in place (SEC, RCI, Multistar, RGLM, Comb)
5. [x] s-t path support alongside closed tours
6. [x] Project/package rename to `rcspp` and pip-installable Python package
7. [x] Sub-MIP SEC separation with column mapping
8. [x] Deterministic warm-start/primal heuristic mode
9. [x] `thread_local` callback handling for parallel solve safety
10. [x] Pluggable `SeparationOracle` and `BoundPropagator` with optional HiGHS build
11. [x] LP-guided primal heuristic callback during solve
12. [x] Shortest Path Inequality (SPI) separator with lifting
13. [x] Reduced-cost fixing and prove-only mode
14. [x] Documentation overhaul and published docs site
15. [ ] Remote branch cleanup for merged/obsolete work
16. [ ] Parameter tuning campaign (separation, propagation, heuristic defaults)
17. [ ] Close the remaining hard Roberti M-series instances

---

## 3. Implementation Roadmap

Each **work unit** (e.g., 2.1, 3.3, 5.5) is one PR. When you say "do 3.3",
that means: create branch `3.3-<short-name>`, implement it, open a PR.

Work units have explicit file ownership and dependencies. Agents check open
branches and PRs before claiming a unit.

### Step 0 — ng-paths / DSSR / bucket-inspired bounds

```
Deliverable: stronger RCSPP bounds pipeline using iterative ng growth, faster
dominance, and asynchronous bound updates during B&B.
```

| ID | Status | PR title | Deliverable | Files | Depends on |
|----|--------|----------|-------------|-------|------------|
| 0.1 | [x] | Cross-project design sync (baldes + pathwyse) | Import dominance/data-layout ideas from `../blades` and DSSR growth strategy from `../pathwyse`; document target architecture and invariants | docs/, src/preprocess/ | — |
| 0.2 | [x] | SoA + SIMD dominance baseline | Refactor labeling dominance path to SoA layout with AVX/scalar fallback; add deterministic microbench and correctness parity tests | src/preprocess/, tests/ | 0.1 |
| 0.3 | [x] | Shared bound store for async updates | Add versioned shared store for forward/backward cost bounds so background DSSR workers can publish snapshots while B&B consumes them safely | src/model/, src/preprocess/ | 0.1 |
| 0.4 | [x] | Parallel DSSR schedule during B&B | Start DSSR with max ng size 4, grow toward 8-12 on background threads, and wire snapshot consumption into propagation/fixing callbacks | src/model/, src/preprocess/, src/util/ | 0.2, 0.3 |
| 0.5 | [x] | Reduced-cost fixing cost/benefit study | Quantify DSSR vs 2-cycle elimination for reduced-cost fixing (fixings gained, runtime, node-count impact), then choose dynamic policy | benchmarks/, docs/benchmarks.md, src/model/ | 0.4 |
| 0.6 | [x] | Incumbent/optimal handoff from ng solver | When ng/DSSR finds full feasible path, inject incumbent; if global proof condition is met, stop B&B and report optimal | src/model/highs_bridge.cpp, src/preprocess/ | 0.4 |
| 0.7 | [x] | Iterative LB reuse and monotone pruning | Reuse DSSR bounds across iterations, enforce monotone tightening, and apply updated bounds in subsequent pruning passes | src/preprocess/, src/model/ | 0.4 |

### Step 1 — Immediate Next Steps

```
Deliverable: deterministic profiling, tuning, and benchmarking baseline.
All runs in this step must be deterministic (fixed seed, fixed thread count,
stable ordering, reproducible machine config).
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 1.1 | Performance analysis and acceleration feasibility | Profile heuristic, separators, and propagator time/memory; assess SIMD/AVX and GPU offload candidates with expected ROI | benchmarks/, docs/ | — |
| 1.2 | Parameter sweep for default settings | Sweep core parameters on medium/large instances that solve within 100s; propose deterministic default values | src/model/, src/heuristic/, src/preprocess/, benchmarks/ | 1.1 |
| 1.3 | Full benchmark + PathWyse comparison | Run complete SPPRCLIB + Roberti benchmark campaign; publish runtime/gap tables and PathWyse comparison | benchmarks/, docs/benchmarks.md | 1.2 |

### Step 2 — Land open PRs

```
Deliverable: merge stalled PRs, clean up branches.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 2.1 | Land deterministic solver (#17) | Review, fix conflicts, merge | src/heuristic/, src/sep/ | — |
| 2.2 | Land thread_local callbacks (#16) | Review, fix conflicts, merge | third_party/highs_patch/ | — |
| 2.3 | Land pluggable separators (#15) | Review, fix conflicts, merge | src/sep/separation_oracle.*, src/preprocess/bound_propagator.* | — |
| 2.4 | Land benchmarks reorg (#13) | Review, fix conflicts, merge | benchmarks/, tests/ | — |
| 2.5 | Land SPI separator | Create PR from branch, review, merge | src/sep/spi_separator.* | — |
| 2.6 | Land cutoff/prove-only | Create PR from branch, review, merge | src/model/ | — |
| 2.7 | Land docs fixes | Create PR from branch, review, merge | docs/ | — |
| 2.8 | Clean up stale branches | Delete merged/abandoned remote branches | — | 2.1–2.7 |

**All of 2.1–2.7 are parallel.** 2.8 depends on all others.

### Step 3 — Test coverage & stability

```
Deliverable: comprehensive test suite, confidence in correctness.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 3.1 | Path vs tour formulation tests | Side-by-side tests: same graph solved as tour and s-t path, verify degree constraints and SEC form | tests/test_formulation.cpp | 2.1 |
| 3.2 | Separator edge case tests | Empty cuts, singleton sets, numerical stability (tight tolerances), all 5 separators | tests/test_separators.cpp | 2.1 |
| 3.3 | Sub-MIP column mapping tests | Verify SEC separation works correctly in sub-MIP context with column mapping | tests/test_separators.cpp | 2.3 |
| 3.4 | Heuristic unit tests | Direct tests for primal heuristic construction + local search (not just via model) | tests/test_heuristic.cpp | 2.1 |
| 3.5 | Python binding integration tests | End-to-end: load instance, solve, check solution from Python | tests/python/test_bindings.py | 2.3 |
| 3.6 | Large instance scaling tests | 200+ node instances, verify solver doesn't degrade unexpectedly | tests/test_instances.cpp | 2.4 |
| 3.7 | Domain propagator regression tests | Edge cases: all-pairs mode, tightened UBs, zero-demand nodes | tests/test_propagator.cpp | 2.1 |

**All of 3.1–3.7 are parallel** — they touch different test files.

### Step 4 — Parameter tuning

```
Deliverable: optimal default parameters for each solver component, backed by
benchmark data. Each work unit profiles one component across the full instance
set and proposes tuned defaults.
All runs must be deterministic (fixed seed/thread count and stable ordering).
```

**Cut separation parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 4.1 | Tune separation_tol | Profile violation threshold (0.01–0.5) across SPPRCLIB + Roberti; find optimal default (current: 0.1, Jepsen et al. suggest 0.4) | src/model/highs_bridge.cpp, benchmarks/ | 2.4 |
| 4.2 | Tune max_cuts_per_separator | Profile top-k cuts per round (1–20) per separator; may differ per separator family | src/model/highs_bridge.cpp, benchmarks/ | 2.4 |
| 4.3 | Tune separation_interval | Profile amortization (1–5) across instance size classes; larger instances may benefit from less frequent separation | src/model/highs_bridge.cpp, benchmarks/ | 2.4 |
| 4.4 | RGLM enable/disable decision | Benchmark RGLM on full instance set with tuned params from 4.1–4.3; decide if it should be on by default | src/sep/rglm_separator.cpp, benchmarks/ | 4.1, 4.2 |
| 4.5 | Comb effectiveness evaluation | Benchmark Comb on full instance set; if ineffective, deprecate with documented rationale | src/sep/comb_separator.cpp, benchmarks/ | 4.1, 4.2 |

**Propagator parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 4.6 | Tune propagator label budget | Profile max labels per node (10–200, current hardcoded 50) on hard instances | src/preprocess/edge_elimination.h, benchmarks/ | 2.4 |
| 4.7 | Tune all_pairs_propagation | Profile all-pairs mode cost vs fixing gain; find instance-size threshold where it pays off | src/model/model.cpp, benchmarks/ | 2.4 |

**Heuristic parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 4.8 | Tune primal heuristic (warm-start) time budget | Profile time budget scaling (current: min(500ms, n*10ms)) vs solution quality and total solve time | src/heuristic/warm_start.h, benchmarks/ | 2.1 |
| 4.9 | Tune primal heuristic (warm-start) restart count | Profile number of restarts (current: time-based) vs quality; ties into deterministic solver (#17) | src/heuristic/warm_start.h, benchmarks/ | 2.1 |
| 4.10 | Tune local search neighbourhood | Profile 2-opt vs or-opt vs combined; measure move acceptance rates and time per move | src/heuristic/warm_start.h, benchmarks/ | 2.1 |

**4.1, 4.2, 4.3 are parallel** (independent parameter sweeps on separation).
**4.6, 4.7 are parallel** (independent propagator experiments).
**4.8, 4.9, 4.10 are parallel** (independent heuristic experiments).
4.4 and 4.5 depend on 4.1 + 4.2. All three groups are parallel with each other.

### Step 5 — Solver internals & performance

```
Deliverable: tackle the 5 unsolved Roberti M-series instances (150–200 nodes).
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 5.1 | CPU/memory profiling on hard instances | Identify bottlenecks: LP vs separation vs propagation | benchmarks/ | 2.4 |
| 5.2 | Branching strategy control | Custom variable selection via HiGHS callbacks (if available) | src/model/highs_bridge.cpp | 2.2 |
| 5.3 | Parallel Model::solve() | Enable concurrent solves on different threads (uses thread_local from 2.2) | src/model/ | 2.2 |
| 5.4 | Early termination heuristic | Skip separation when UB stabilizes or gap is tiny | src/model/highs_bridge.cpp | 2.1 |
| 5.5 | Ryan-Foster branching for RCSPP | Investigate branching on edge pairs (from threading analysis branch) | src/model/ | 5.2 |

**5.1, 5.2, 5.3, 5.4 are parallel.** 5.5 depends on 5.2.

---

### Step 6 — Applications Expansion

```
Deliverable: prioritize use-case support from docs/applications.md using
published papers and public instance sets.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 6.1 | PCTSP / Selective TSP onboarding | Consolidate academic references, ingest/normalize public instances, document objective/penalty convention mapping to current CPTP model | docs/applications.md, docs/benchmarks.md, src/core/ | 1.3 |
| 6.2 | Orienteering edge-resource extension | Add support plan for edge-budget resources and instance ingestion path for OP datasets | docs/applications.md, src/model/, src/core/ | 6.1 |
| 6.3 | Multi-resource extension planning | Define two-resource formulation/test plan and identify benchmark sets for capacity+time style pricing | docs/applications.md, docs/roadmap.md | 6.2 |
| 6.4 | Directed support planning (last) | Define directed-formulation migration path and collect asymmetric benchmark instances | docs/applications.md, docs/roadmap.md | 6.3 |

---
