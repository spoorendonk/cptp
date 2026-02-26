# rcspp-bac — Implementation Roadmap

Branch-and-cut solver for the Resource Constrained Shortest Path Problem (RCSPP).

Current state: production-quality solver with SEC, RCI, Multistar, RGLM, Comb
separators, demand-reachability preprocessing, edge elimination via labeling,
domain propagation, and parallel primal heuristic (warm-start). Solves 45/45 SPPRCLIB
instances and 26/31 Roberti Set 3 instances to optimality.

---

## 1. Current Architecture

```
src/core/        — Problem, static_graph, IO, Dinitz max-flow, Gomory-Hu tree
src/preprocess/  — Demand-reachability, edge elimination (label-correcting)
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

Each **work unit** (e.g., 1.1, 2.3, 4.5) is one PR. When you say "do 2.3",
that means: create branch `2.3-<short-name>`, implement it, open a PR.

Work units have explicit file ownership and dependencies. Agents check open
branches and PRs before claiming a unit.

### Step 0 — Immediate Next Steps

```
Deliverable: deterministic profiling, tuning, and benchmarking baseline.
All runs in this step must be deterministic (fixed seed, fixed thread count,
stable ordering, reproducible machine config).
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 0.1 | Performance analysis and acceleration feasibility | Profile heuristic, separators, and propagator time/memory; assess SIMD/AVX and GPU offload candidates with expected ROI | benchmarks/, docs/ | — |
| 0.2 | Parameter sweep for default settings | Sweep core parameters on medium/large instances that solve within 100s; propose deterministic default values | src/model/, src/heuristic/, src/preprocess/, benchmarks/ | 0.1 |
| 0.3 | Full benchmark + PathWyse comparison | Run complete SPPRCLIB + Roberti benchmark campaign; publish runtime/gap tables and PathWyse comparison | benchmarks/, docs/benchmarks.md | 0.2 |

### Step 1 — Land open PRs

```
Deliverable: merge stalled PRs, clean up branches.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 1.1 | Land deterministic solver (#17) | Review, fix conflicts, merge | src/heuristic/, src/sep/ | — |
| 1.2 | Land thread_local callbacks (#16) | Review, fix conflicts, merge | third_party/highs_patch/ | — |
| 1.3 | Land pluggable separators (#15) | Review, fix conflicts, merge | src/sep/separation_oracle.*, src/preprocess/bound_propagator.* | — |
| 1.4 | Land benchmarks reorg (#13) | Review, fix conflicts, merge | benchmarks/, tests/ | — |
| 1.5 | Land SPI separator | Create PR from branch, review, merge | src/sep/spi_separator.* | — |
| 1.6 | Land cutoff/prove-only | Create PR from branch, review, merge | src/model/ | — |
| 1.7 | Land docs fixes | Create PR from branch, review, merge | docs/ | — |
| 1.8 | Clean up stale branches | Delete merged/abandoned remote branches | — | 1.1–1.7 |

**All of 1.1–1.7 are parallel.** 1.8 depends on all others.

### Step 2 — Test coverage & stability

```
Deliverable: comprehensive test suite, confidence in correctness.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 2.1 | Path vs tour formulation tests | Side-by-side tests: same graph solved as tour and s-t path, verify degree constraints and SEC form | tests/test_formulation.cpp | 1.1 |
| 2.2 | Separator edge case tests | Empty cuts, singleton sets, numerical stability (tight tolerances), all 5 separators | tests/test_separators.cpp | 1.1 |
| 2.3 | Sub-MIP column mapping tests | Verify SEC separation works correctly in sub-MIP context with column mapping | tests/test_separators.cpp | 1.3 |
| 2.4 | Heuristic unit tests | Direct tests for primal heuristic construction + local search (not just via model) | tests/test_heuristic.cpp | 1.1 |
| 2.5 | Python binding integration tests | End-to-end: load instance, solve, check solution from Python | tests/python/test_bindings.py | 1.3 |
| 2.6 | Large instance scaling tests | 200+ node instances, verify solver doesn't degrade unexpectedly | tests/test_instances.cpp | 1.4 |
| 2.7 | Domain propagator regression tests | Edge cases: all-pairs mode, tightened UBs, zero-demand nodes | tests/test_propagator.cpp | 1.1 |

**All of 2.1–2.7 are parallel** — they touch different test files.

### Step 3 — Parameter tuning

```
Deliverable: optimal default parameters for each solver component, backed by
benchmark data. Each work unit profiles one component across the full instance
set and proposes tuned defaults.
All runs must be deterministic (fixed seed/thread count and stable ordering).
```

**Cut separation parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 3.1 | Tune separation_tol | Profile violation threshold (0.01–0.5) across SPPRCLIB + Roberti; find optimal default (current: 0.1, Jepsen et al. suggest 0.4) | src/model/highs_bridge.cpp, benchmarks/ | 1.4 |
| 3.2 | Tune max_cuts_per_separator | Profile top-k cuts per round (1–20) per separator; may differ per separator family | src/model/highs_bridge.cpp, benchmarks/ | 1.4 |
| 3.3 | Tune separation_interval | Profile amortization (1–5) across instance size classes; larger instances may benefit from less frequent separation | src/model/highs_bridge.cpp, benchmarks/ | 1.4 |
| 3.4 | RGLM enable/disable decision | Benchmark RGLM on full instance set with tuned params from 3.1–3.3; decide if it should be on by default | src/sep/rglm_separator.cpp, benchmarks/ | 3.1, 3.2 |
| 3.5 | Comb effectiveness evaluation | Benchmark Comb on full instance set; if ineffective, deprecate with documented rationale | src/sep/comb_separator.cpp, benchmarks/ | 3.1, 3.2 |

**Propagator parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 3.6 | Tune propagator label budget | Profile max labels per node (10–200, current hardcoded 50) on hard instances | src/preprocess/edge_elimination.h, benchmarks/ | 1.4 |
| 3.7 | Tune all_pairs_propagation | Profile all-pairs mode cost vs fixing gain; find instance-size threshold where it pays off | src/model/model.cpp, benchmarks/ | 1.4 |

**Heuristic parameters:**

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 3.8 | Tune primal heuristic (warm-start) time budget | Profile time budget scaling (current: min(500ms, n*10ms)) vs solution quality and total solve time | src/heuristic/warm_start.h, benchmarks/ | 1.1 |
| 3.9 | Tune primal heuristic (warm-start) restart count | Profile number of restarts (current: time-based) vs quality; ties into deterministic solver (#17) | src/heuristic/warm_start.h, benchmarks/ | 1.1 |
| 3.10 | Tune local search neighbourhood | Profile 2-opt vs or-opt vs combined; measure move acceptance rates and time per move | src/heuristic/warm_start.h, benchmarks/ | 1.1 |

**3.1, 3.2, 3.3 are parallel** (independent parameter sweeps on separation).
**3.6, 3.7 are parallel** (independent propagator experiments).
**3.8, 3.9, 3.10 are parallel** (independent heuristic experiments).
3.4 and 3.5 depend on 3.1 + 3.2. All three groups are parallel with each other.

### Step 4 — Solver internals & performance

```
Deliverable: tackle the 5 unsolved Roberti M-series instances (150–200 nodes).
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 4.1 | CPU/memory profiling on hard instances | Identify bottlenecks: LP vs separation vs propagation | benchmarks/ | 1.4 |
| 4.2 | Branching strategy control | Custom variable selection via HiGHS callbacks (if available) | src/model/highs_bridge.cpp | 1.2 |
| 4.3 | Parallel Model::solve() | Enable concurrent solves on different threads (uses thread_local from 1.2) | src/model/ | 1.2 |
| 4.4 | Early termination heuristic | Skip separation when UB stabilizes or gap is tiny | src/model/highs_bridge.cpp | 1.1 |
| 4.5 | Ryan-Foster branching for RCSPP | Investigate branching on edge pairs (from threading analysis branch) | src/model/ | 4.2 |

**4.1, 4.2, 4.3, 4.4 are parallel.** 4.5 depends on 4.2.

---

### Step 5 — Applications Expansion

```
Deliverable: prioritize use-case support from docs/applications.md using
published papers and public instance sets.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 5.1 | PCTSP / Selective TSP onboarding | Consolidate academic references, ingest/normalize public instances, document objective/penalty convention mapping to current CPTP model | docs/applications.md, docs/benchmarks.md, src/core/ | 0.3 |
| 5.2 | Orienteering edge-resource extension | Add support plan for edge-budget resources and instance ingestion path for OP datasets | docs/applications.md, src/model/, src/core/ | 5.1 |
| 5.3 | Multi-resource extension planning | Define two-resource formulation/test plan and identify benchmark sets for capacity+time style pricing | docs/applications.md, docs/roadmap.md | 5.2 |
| 5.4 | Directed support planning (last) | Define directed-formulation migration path and collect asymmetric benchmark instances | docs/applications.md, docs/roadmap.md | 5.3 |

---
