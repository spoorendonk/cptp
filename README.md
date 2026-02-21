# rcspp-bac

Branch-and-cut solver for a variant of the **Resource Constrained Shortest Path Problem (RCSPP)**.

## Problem Variant

The solver targets an RCSPP with:

- **Single capacity resource**: each vertex has a demand; the total demand on the path must not exceed vehicle capacity Q
- **Vertex profits**: each vertex has a profit collected when visited; the objective minimizes `travel_cost - collected_profit`
- **Optional vertices**: not all vertices need to be visited -- the solver selects the optimal subset
- **Negative-cost cycles**: the graph may contain negative cycles (profits can exceed edge costs), requiring elementary path constraints
- **Tours or s-t paths**: supports closed tours (depot-to-depot) and open s-t paths (source-to-target) with automatic detection

This problem is also known as the **Capacitated Profitable Tour Problem (CPTP)**, introduced by [Jepsen, Petersen, Spoorendonk & Pisinger (2014)](https://doi.org/10.1016/S1572-5286(14)00036-X). The s-t path variant extends it to open paths where source and target differ.

## Features

- **MIP formulation** with binary edge and node variables (undirected)
- **Dynamic cut separation**: subtour elimination (SEC), rounded capacity inequalities (RCI), multistar inequalities
- **Gomory-Hu tree** (Gusfield's algorithm) shared across separators for efficient min-cut computation
- **Preprocessing**: demand-reachability filtering and ESPPRC-style edge elimination
- **Warm-start heuristic**: parallel greedy construction + local search (2-opt, or-opt, node drop/add) via TBB
- **MIP backend**: HiGHS with user cut callbacks
- **Python bindings** via nanobind (optional)

## Build

Requires GCC 14+ (C++23), CMake 3.25+, and TBB:

```bash
apt install libtbb-dev
git submodule update --init
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

For Python bindings:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCPTP_BUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

## Usage

```bash
./build/cptp-solve <instance> [--source <node>] [--target <node>] [--time_limit <sec>] [--threads <n>]
```

Accepts TSPLIB (`.vrp`, `.sppcc`) and PathWyse (`.txt`) instance formats. All additional options are forwarded to HiGHS.

When `source != target`, the solver uses an open s-t path formulation (degree 1 at source/target, degree 2 at intermediates). When `source == target` (default), the standard tour formulation is used.

### Examples

```bash
# Tour (closed loop from depot)
./build/cptp-solve bench/instances/spprclib/B-n45-k6-54.sppcc --time_limit 120

# s-t path (open path from node 0 to node 3)
./build/cptp-solve tests/data/tiny4_path.txt

# Override source/target via CLI
./build/cptp-solve tests/data/tiny4.txt --source 0 --target 3
```

## Tests

```bash
./build/cptp_tests
```

## Project Structure

```
src/core/        Problem definition, IO (TSPLIB & PathWyse), Dinitz max-flow, Gomory-Hu tree
src/sep/         Cut separators (SEC, RCI, Multistar)
src/model/       HiGHS MIP integration and cut callbacks
src/heuristic/   Warm-start construction + local search
src/preprocess/  Demand-reachability and edge elimination
src/cli/         Command-line solver
python/          nanobind Python bindings
tests/           Catch2 unit and regression tests
docs/            Algorithm documentation
```

## Documentation

- [Algorithms and techniques](docs/algorithms.md) -- formulation (tours and s-t paths), separators, preprocessing, references
- [Warm-start heuristic](docs/warm-start-heuristic.md) -- construction, local search, parallelism

## Dependencies

- [HiGHS](https://highs.dev) v1.10.0 -- MIP solver (fetched via CMake)
- [melon](https://github.com/fhamonic/melon) -- graph data structures (git submodule)
- [Catch2](https://github.com/catchorg/Catch2) v3.7.1 -- testing (fetched via CMake)
- [TBB](https://github.com/oneapi-src/oneTBB) 2021.11 -- parallel heuristic

## Benchmarks

Tested against [SPPRCLIB](https://or.rwth-aachen.de/research/spprclib) instances and verified against [PathWyse](https://doi.org/10.1080/10556788.2023.2296978).

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). [A branch-and-cut algorithm for the capacitated profitable tour problem](https://doi.org/10.1016/S1572-5286(14)00036-X). *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
- Salani, M., Basso, S., & Giuffrida, V. (2024). [PathWyse: a flexible, open-source library for the resource constrained shortest path problem](https://doi.org/10.1080/10556788.2023.2296978). *Optimization Methods and Software*, 39(2).
