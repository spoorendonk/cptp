# rcspp-bac — Resource Constrained Shortest Path Branch-and-Cut Solver

## Git Workflow

- **Never commit directly to main.** Always create a feature branch, push, and open a PR.
- **Linear history only.** Merge PRs with squash or rebase (no merge commits).
- **No force-push to main.**

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
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
./build/rcspp_tests                              # C++ tests
./build/rcspp-solve tests/data/tiny4.txt         # CLI
```

## Dependencies

- GCC 14, C++23
- CMake 3.25+
- `apt install libtbb-dev` (TBB 2021.11)
- HiGHS: fetched via CMake FetchContent (v1.10.0)
- Catch2: fetched via CMake FetchContent (v3.7.1)

## Architecture

```
src/core/        — Problem definition, IO parsers (TSPLIB, numeric .txt), Dinitz max-flow, Gomory-Hu tree
src/preprocess/  — Demand-reachability and edge elimination via capacity-aware labeling
src/sep/         — Solver-independent separators (SEC, RCI, Multistar, RGLM, Comb)
src/model/       — HiGHS integration (Model, HiGHSBridge, propagator)
src/heuristic/   — Primal heuristics (initial solution + LP-guided callback)
src/cli/         — CLI tool (rcspp-solve)
src/util/        — Utilities (Timer)
python/          — nanobind Python bindings
tests/           — Catch2 tests
docs/            — Algorithm documentation
```

## Key types

- `rcspp::Problem` — RCSPP instance using `static_graph` (own CSR graph); has `source()`, `target()`, `is_tour()`
- `rcspp::Model` — User-facing solver interface; `set_source()`/`set_target()` for paths, `set_depot()` for tours
- `rcspp::HiGHSBridge` — Wires separators + domain propagator into HiGHS MIP
- `rcspp::sep::Separator` — Base class for cut separators
- `rcspp::sep::SECSeparator` — Subtour elimination via Dinitz max-flow (path-aware)
- `rcspp::gomory_hu_tree` — Gusfield's algorithm, shared across separators
- `rcspp::heuristic::build_initial_solution` — Pre-solve greedy + local search heuristic
- `rcspp::heuristic::lp_guided_heuristic` — LP-guided callback heuristic (reduced graph)

## Tour vs s-t path

When `source == target` (default), the solver uses a closed tour formulation (degree 2 at all nodes).
When `source != target`, it uses an open path formulation:
- Degree 1 at source/target, degree 2 at intermediates
- SEC cuts: sets containing the path target need only 1 cut crossing (path enters and terminates)
- Numeric `.txt` format: optional `source target` line after the capacity line

## Namespace

All code under `rcspp::` namespace, separators under `rcspp::sep::`, heuristics under `rcspp::heuristic::`, preprocessing under `rcspp::preprocess::`.
