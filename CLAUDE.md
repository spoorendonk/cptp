# rcspp-bac — Resource Constrained Shortest Path Branch-and-Cut Solver

## Git Workflow

- **Never commit directly to main.** Always create a feature branch, push, and open a PR.
- **If user says "commit" while on main**: create a feature branch, commit there, push, and open a PR automatically.
- **Linear history only.** Merge PRs with squash or rebase (no merge commits).
- **No force-push to main.**

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Algorithms only (no HiGHS):
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_HIGHS=OFF
cmake --build build -j$(nproc)
```

For Python bindings (development):
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

For pip-installable Python package:
```bash
pip install .
```

## Test

```bash
./build/rcspp_algo_tests                         # 57 algorithm tests (no HiGHS needed)
./build/rcspp_tests                              # 26 integration tests (requires HiGHS)
pytest tests/python/test_algorithms.py           # 33 Python algorithm tests
pytest tests/python/test_solver.py               # Python solver tests (requires HiGHS)
./build/rcspp-solve tests/data/tiny4.txt         # CLI
```

Test structure:
- `tests/test_max_flow.cpp` — Dinitz max-flow, Gomory-Hu tree (5 tests)
- `tests/test_separators.cpp` — Individual separators (SEC, RCI, Multistar, Comb, RGLM), preprocessing, warm-start heuristic (29 tests)
- `tests/test_oracle.cpp` — SeparationOracle, BoundPropagator pluggable API (23 tests)
- `tests/test_model.cpp` — Model API, solver integration [requires HiGHS]
- `tests/test_optimal.cpp` — Known-optimal instances from SPPRCLIB [requires HiGHS]
- `tests/test_propagator.cpp` — Labeling, edge elimination, propagator integration [requires HiGHS]
- `tests/python/test_algorithms.py` — Python bindings for algorithms (33 tests, no HiGHS)
- `tests/python/test_solver.py` — Python solver API [requires HiGHS]

## Dependencies

- GCC 14, C++23
- CMake 3.25+
- `apt install libtbb-dev` (TBB 2021.11)
- HiGHS: fetched via CMake FetchContent (v1.10.0)
- Catch2: fetched via CMake FetchContent (v3.7.1)

## Architecture

Two libraries:
- **`rcspp_algorithms`** — solver-independent (separators, propagation, heuristic, preprocessing). Only TBB.
- **`rcspp_model`** — HiGHS integration (optional, `RCSPP_BUILD_HIGHS=ON` by default).

```
src/core/        — Problem definition, IO parsers (TSPLIB, numeric .txt), Dinitz max-flow, Gomory-Hu tree
src/sep/         — Solver-independent separators (SEC, RCI, Multistar, RGLM, Comb) + SeparationOracle
src/preprocess/  — BoundPropagator, demand-reachability, edge elimination via capacity-aware labeling
src/heuristic/   — Warm-start construction + local search
src/model/       — HiGHS integration (Model, HiGHSBridge, propagator) [optional]
src/cli/         — CLI tool (rcspp-solve) [requires HiGHS]
src/util/        — Utilities (Logger, Timer)
python/          — nanobind Python bindings
tests/           — Catch2 tests + Python tests
benchmarks/      — Benchmark instances, scripts, and results
docs/            — Algorithm documentation
```

## Key types

- `rcspp::Problem` — RCSPP instance using `static_graph` (own CSR graph); has `source()`, `target()`, `is_tour()`
- `rcspp::sep::SeparationOracle` — Solver-independent cut separation (bundles support graph, GH tree, parallel separators)
- `rcspp::preprocess::BoundPropagator` — Solver-independent domain propagation (sweep + chain fixings)
- `rcspp::heuristic::build_warm_start` — Parallel greedy + local search heuristic
- `rcspp::sep::Separator` — Base class for cut separators
- `rcspp::sep::SECSeparator` — Subtour elimination via Dinitz max-flow (path-aware)
- `rcspp::gomory_hu_tree` — Gusfield's algorithm, shared across separators
- `rcspp::Model` — User-facing solver interface (HiGHS); `set_source()`/`set_target()` for paths, `set_depot()` for tours
- `rcspp::HiGHSBridge` — Wires separators + domain propagator into HiGHS MIP

## Tour vs s-t path

When `source == target` (default), the solver uses a closed tour formulation (degree 2 at all nodes).
When `source != target`, it uses an open path formulation:
- Degree 1 at source/target, degree 2 at intermediates
- SEC cuts: sets containing the path target need only 1 cut crossing (path enters and terminates)
- Numeric `.txt` format: optional `source target` line after the capacity line

## Namespace

All code under `rcspp::` namespace, separators under `rcspp::sep::`, heuristics under `rcspp::heuristic::`, preprocessing under `rcspp::preprocess::`.

## Agent Coordination

The roadmap (docs/ROADMAP.md) defines work units with IDs like `2.3`,
`3.8`, `5.1`. Each work unit maps to a branch name (e.g., `3.1-tune-separation-tol`).

**Before suggesting or starting any work unit:**
1. Check open branches: `git branch -a`
2. Check open PRs: `gh pr list`
3. Check for running agents on this machine (background tasks, worktrees)
4. Never start a work unit that another agent has an open branch or PR for
5. Prefer the lowest-numbered unblocked, unclaimed work unit

## Fullgate

When the user says **"fullgate"**, run this sequence in order. Each step can also be invoked individually by name:

1. **Feature branch** — create one if not already on a feature branch
2. **Create PR** — if no PR exists for the current branch
3. **Sync main** — pull latest main and merge into the current feature branch, resolve conflicts
4. **Tests** — check if new/updated tests are needed and add them
5. **Update docs** — update README.md and docs as needed
6. **Push & update PR**
7. **Review** — thoroughly review the PR (code quality, correctness, style, tests, performance)
8. **Build** — `cmake --build build -j$(nproc)`
9. **Test** — `./build/rcspp_tests`
10. **Push & update PR** again with any fixes
11. **Finalize** — if nothing more to do: squash-merge the PR, delete feature branch (local + remote), pull main, switch to main
