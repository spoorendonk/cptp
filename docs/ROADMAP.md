# rcspp-bac — Implementation Roadmap

Branch-and-cut solver for the Resource Constrained Shortest Path Problem (RCSPP).

Current state: production-quality solver with SEC, RCI, Multistar, RGLM, Comb
separators, demand-reachability preprocessing, edge elimination via labeling,
domain propagation, and parallel warm-start heuristic. Solves 45/45 SPPRCLIB
instances and 26/31 Roberti Set 3 instances to optimality.

---

## 1. Current Architecture

```
src/core/        — Problem, static_graph, IO, Dinitz max-flow, Gomory-Hu tree
src/preprocess/  — Demand-reachability, edge elimination (label-correcting)
src/sep/         — SEC, RCI, Multistar, RGLM, Comb separators
src/model/       — HiGHS integration (Model, HiGHSBridge, propagator)
src/heuristic/   — Warm-start construction + local search
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

## 2. Open Work (from branches and PRs)

These branches/PRs represent in-progress or stalled work. Check status before
starting any overlapping work unit.

| Branch | PR | Status | Summary |
|--------|-----|--------|---------|
| `claude/make-deterministic-ZwXew` | #17 | Open | Deterministic warm-start (fixed restarts, seeded RNG, stable_sort) |
| `claude/extract-highs-patching-kWX1I` | #16 | Open | thread_local callbacks for parallel Model::solve() |
| `claude/pluggable-separators-interface-nRwPT` | #15 | Open | SeparationOracle + BoundPropagator, optional HiGHS build |
| `reorganize-benchmarks` | #13 | Done | Rename bench/ → benchmarks/, PathWyse workflow |
| `claude/shortest-path-inequalities-ID51D` | — | Branch | Shortest path inequality (SPI) separator |
| `claude/test-ubs-cutoff-performance-FWIOd` | — | Branch | Cutoff + prove-only mode options |
| `claude/github-pages-exploration-7DNje` | — | Branch | MkDocs Material documentation site |
| `claude/highs-threading-analysis-yKBB5` | — | Branch | Ryan-Foster branching, parallel MIP design |
| `claude/review-docs-links-75ohk` | — | Branch | Doc fixes (Roberti attribution, LaTeX math) |

---

## 3. Implementation Roadmap

Each **work unit** (e.g., 1.1, 2.3, 5.7) is one PR. When you say "do 2.3",
that means: create branch `2.3-<short-name>`, implement it, open a PR.

Work units have explicit file ownership and dependencies. Agents check open
branches and PRs before claiming a unit.

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
| 2.4 | Heuristic unit tests | Direct tests for warm-start construction + local search (not just via model) | tests/test_heuristic.cpp | 1.1 |
| 2.5 | Python binding integration tests | End-to-end: load instance, solve, check solution from Python | tests/python/test_bindings.py | 1.3 |
| 2.6 | Large instance scaling tests | 200+ node instances, verify solver doesn't degrade unexpectedly | tests/test_instances.cpp | 1.4 |
| 2.7 | Domain propagator regression tests | Edge cases: all-pairs mode, tightened UBs, zero-demand nodes | tests/test_propagator.cpp | 1.1 |

**All of 2.1–2.7 are parallel** — they touch different test files.

### Step 3 — Parameter tuning

```
Deliverable: optimal default parameters for each solver component, backed by
benchmark data. Each work unit profiles one component across the full instance
set and proposes tuned defaults.
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
| 3.8 | Tune warm-start time budget | Profile time budget scaling (current: min(500ms, n*10ms)) vs solution quality and total solve time | src/heuristic/warm_start.h, benchmarks/ | 1.1 |
| 3.9 | Tune warm-start restart count | Profile number of restarts (current: time-based) vs quality; ties into deterministic solver (#17) | src/heuristic/warm_start.h, benchmarks/ | 1.1 |
| 3.10 | Tune local search neighbourhood | Profile 2-opt vs or-opt vs combined; measure move acceptance rates and time per move | src/heuristic/warm_start.h, benchmarks/ | 1.1 |

**3.1, 3.2, 3.3 are parallel** (independent parameter sweeps on separation).
**3.6, 3.7 are parallel** (independent propagator experiments).
**3.8, 3.9, 3.10 are parallel** (independent heuristic experiments).
3.4 and 3.5 depend on 3.1 + 3.2. All three groups are parallel with each other.

### Step 4 — Python builds & CI

```
Deliverable: CI pipeline, cross-platform wheels, PyPI-ready releases.
Follows mip-heuristics-2 pattern: scikit-build-core + nanobind + cibuildwheel.
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 4.1 | CI workflow for C++ build + test | GitHub Actions: build with g++-14, run rcspp_tests on ubuntu-24.04 | .github/workflows/ci.yml | 1.1 |
| 4.2 | CI workflow for Python build + test | GitHub Actions: `pip install -v .`, run pytest on ubuntu-24.04 | .github/workflows/ci.yml, tests/python/ | 2.5 |
| 4.3 | Add _version.py and stable ABI | Single version source in python/rcspp_bac/_version.py, add `wheel.py-api = "cp312"` to pyproject.toml | python/rcspp_bac/_version.py, pyproject.toml | 1.1 |
| 4.4 | cibuildwheel config | Build wheels for Linux x86_64/aarch64, macOS arm64, Windows amd64 (cp312+) | pyproject.toml | 4.3 |
| 4.5 | Release workflow | GitHub Actions on tag push: cibuildwheel → sdist → publish to PyPI | .github/workflows/release.yml | 4.4 |

**4.1, 4.3 are parallel.** 4.2 depends on 2.5. 4.4 depends on 4.3. 4.5 depends on 4.4.

### Step 5 — Solver internals & performance

```
Deliverable: tackle the 5 unsolved Roberti M-series instances (150–200 nodes).
```

| ID | PR title | Deliverable | Files | Depends on |
|----|----------|-------------|-------|------------|
| 5.1 | CPU/memory profiling on hard instances | Identify bottlenecks: LP vs separation vs propagation | benchmarks/ | 1.4 |
| 5.2 | Branching strategy control | Custom variable selection via HiGHS callbacks (if available) | src/model/highs_bridge.cpp | 1.2 |
| 5.3 | Parallel Model::solve() | Enable concurrent solves on different threads (uses thread_local from 1.2) | src/model/ | 1.2 |
| 5.4 | Early termination heuristic | Skip separation when UB stabilizes or gap is tiny | src/model/highs_bridge.cpp | 1.1 |
| 5.5 | Ryan-Foster branching for RCSPP | Investigate branching on edge pairs (from threading analysis branch) | src/model/ | 5.2 |

**5.1, 5.2, 5.3, 5.4 are parallel.** 5.5 depends on 5.2.

---

## 4. Parallelism Summary

```
Step 1: 7 parallel work units (all independent PR reviews)
Step 2: 7 parallel work units (all different test files)
Step 3: 3 groups of parallel sweeps (cuts, propagator, heuristic) — 10 total
Step 4: 2 parallel lanes (CI + version), then wheels, then release
Step 5: 4 parallel, then Ryan-Foster
```

**Cross-step parallelism:**
- Steps 2 and 3 are fully parallel (tests vs parameter tuning)
- Step 4 overlaps with step 3 (CI setup while tuning runs)
- Step 5 is parallel with steps 3 and 4 (solver internals vs tuning vs CI)

---

## 5. Agent Coordination

Each work unit ID (e.g., "3.1") maps to a branch name (`3.1-tune-separation-tol`).

Before starting a work unit:
1. `git branch -a` — check for open branches
2. `gh pr list` — check for open PRs
3. Never start a work unit that another agent has an open branch for
4. Prefer the lowest-numbered unblocked, unclaimed work unit

---

## 6. Priority Assessment

**Highest ROI (do first):**
- Step 1 (land open PRs) — unlocks everything else
- Step 2 (tests) — confidence before changes
- Step 3 (parameter tuning) — data-driven defaults, low risk, high impact

**Medium ROI:**
- Step 4 (Python builds) — CI + release pipeline
- Step 5.1 (profiling) — identify where time is spent

**Long-term:**
- Step 5.2–5.5 (solver internals) — guided by profiling results

**Frontier target:** solve the 5 remaining Roberti M-series instances (151–200 nodes).
Steps 3 and 5 contribute directly to this goal.
