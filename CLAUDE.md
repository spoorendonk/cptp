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

For Python bindings:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCPTP_BUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

## Test

```bash
./build/cptp_tests                              # C++ tests
./build/cptp-solve tests/data/tiny4.txt         # CLI
```

## Dependencies

- GCC 14, C++23
- CMake 3.25+
- `apt install libtbb-dev` (TBB 2021.11)
- melon: git submodule in `third_party/melon`
- HiGHS: fetched via CMake FetchContent (v1.10.0)
- Catch2: fetched via CMake FetchContent (v3.7.1)

## Architecture

```
src/core/        — Problem definition, IO parsers (TSPLIB, PathWyse), Dinitz max-flow, Gomory-Hu tree
src/sep/         — Solver-independent separators (SEC, RCI, Multistar, RGLM, Comb)
src/model/       — HiGHS integration (Model, HiGHSBridge)
src/heuristic/   — Warm-start construction + local search
src/preprocess/  — Demand-reachability and edge elimination
src/cli/         — CLI tool (cptp-solve)
src/util/        — Utilities (Timer)
python/          — nanobind Python bindings
tests/           — Catch2 tests
docs/            — Algorithm documentation
```

## Key types

- `cptp::Problem` — CPTP instance using `melon::static_digraph`
- `cptp::Model` — User-facing solver interface
- `cptp::HiGHSBridge` — Wires separators into HiGHS MIP
- `cptp::sep::Separator` — Base class for cut separators
- `cptp::sep::SECSeparator` — Subtour elimination via Dinitz max-flow
- `cptp::gomory_hu_tree` — Gusfield's algorithm, shared across separators
- `cptp::heuristic::build_warm_start` — Parallel greedy + local search heuristic

## Namespace

All code under `cptp::` namespace, separators under `cptp::sep::`, heuristics under `cptp::heuristic::`, preprocessing under `cptp::preprocess::`.
