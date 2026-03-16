# cptp — Capacitated Profitable Tour Problem Solver

## Quick Reference

```bash
# build
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)

# debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build -j$(nproc)

# python editable install
.venv/bin/pip install -e .

# format
clang-format -i src/**/*.{h,cpp}

# test
./build/cptp_tests && .venv/bin/python -m pytest tests/python/test_solver.py

# all tests (including extra: max flow, oracle, separators)
./build/cptp_tests && ./build/cptp_tests_extra && .venv/bin/python -m pytest tests/python/test_solver.py

# benchmarks
benchmarks/run_benchmarks.sh
```

The build patches HiGHS with custom user-cut and propagator callbacks from `third_party/highs_patch/`.

## Git

- Never commit directly to `main`. Always feature branches.
- Linear history only (squash-merge or rebase-merge).
