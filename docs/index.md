---
hide:
  - navigation
---

# rcspp-bac

`rcspp-bac` is a branch-and-cut solver for RCSPP/CPTP variants using HiGHS.

## Getting Started (CLI)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/rcspp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
```

Continue with:
- [CLI Usage](cli.md)
- [Instance Formats](instance-formats.md)
- [API (C++ and Python)](api.md)

## Technical Docs

- [Algorithms Overview](algorithms.md)
- [Cut Separation](separation.md)
- [Preprocessing](preprocessing.md)
- [Domain Propagator](domain-propagator.md)
- [Hyperplane Branching](hyperplane-branching.md)
- [Solve Concurrency](solve-concurrency.md)
- [Primal Heuristic](primal-heuristic.md)
- [ng/DSSR Bounds](ng-dssr-labelling.md)
