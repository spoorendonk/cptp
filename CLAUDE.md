# CPTP Branch-and-Cut Solver

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
src/core/     — Problem definition, IO parsers (TSPLIB, PathWyse)
src/sep/      — Solver-independent separators (SEC, RCI, Multistar, Comb)
src/model/    — HiGHS integration (Model, HiGHSBridge)
src/cli/      — CLI tool (cptp-solve)
src/util/     — Utilities (Timer)
python/       — nanobind Python bindings
tests/        — Catch2 tests
```

## Key types

- `cptp::Problem` — CPTP instance using `melon::static_digraph`
- `cptp::Model` — User-facing solver interface
- `cptp::HiGHSBridge` — Wires separators into HiGHS MIP
- `cptp::sep::Separator` — Base class for cut separators
- `cptp::sep::SECSeparator` — Subtour elimination via Dinitz max-flow

## Namespace

All code under `cptp::` namespace, separators under `cptp::sep::`.
