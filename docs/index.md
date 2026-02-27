---
hide:
  - navigation
---

# rcspp-bac

`rcspp-bac` is a branch-and-cut solver for the Resource Constrained Shortest Path Problem (RCSPP).

## Public Interface

- C++: `rcspp::Model`, `rcspp::Problem`, `rcspp::SolveResult`
- Python: `Model`, `solve`, `Problem`, `SolveResult`, `Status`, `load`
- CLI: `rcspp-solve`

The repository keeps separation, preprocessing, propagation, and heuristic components internal to the solver.

## Quick Start

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/rcspp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
```

## Where to go next

- [Getting Started](getting-started.md)
- [C++ API](cpp-api.md)
- [Python API](python-api.md)
- [Instance Formats](instance-formats.md)
- [Algorithms Overview](algorithms.md)
- [Cut Separation](separation.md)
- [Domain Propagator](domain-propagator.md)
- [Primal Heuristic](primal-heuristic.md)
- [Preprocessing](preprocessing.md)
