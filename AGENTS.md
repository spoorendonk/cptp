# rcspp-bac — Resource Constrained Shortest Path Branch-and-Cut Solver

## Quick Reference

```bash
# build
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)

# build without HiGHS (algorithms only)
cmake -B build -DCMAKE_BUILD_TYPE=Release -DRCSPP_BUILD_HIGHS=OFF && cmake --build build -j$(nproc)

# test
./build/rcspp_algo_tests            # 63 algorithm tests (no HiGHS needed)
./build/rcspp_tests                  # 39 integration tests (requires HiGHS)
pytest tests/python/test_algorithms.py  # 33 Python tests

# python package
pip install .
```

## Git Workflow

- Never commit directly to `main`. Always feature branches.
- If on main when committing: create branch, commit, push, open PR.
- Linear history (squash-merge or rebase-merge). No force-push to `main`.

## Architecture

Two libraries:
- **`rcspp_algorithms`** — solver-independent (separators, propagation, heuristic, preprocessing). Only TBB.
- **`rcspp_model`** — HiGHS integration (optional, default ON).

```
src/core/        — Problem, IO parsers (TSPLIB, numeric .txt), Dinitz max-flow, Gomory-Hu tree
src/sep/         — Separators (SEC, RCI, Multistar, RGLM, Comb) + SeparationOracle
src/preprocess/  — BoundPropagator, demand-reachability, edge elimination
src/heuristic/   — Primal heuristics (initial solution + LP-guided callback)
src/model/       — HiGHS integration (Model, HiGHSBridge, propagator)
src/cli/         — CLI tool (rcspp-solve)
python/          — nanobind Python bindings
tests/           — Catch2 + Python tests
docs/            — Algorithm documentation
```

See `docs/` for algorithm details and theory.

## Key Types

- `rcspp::Problem` — RCSPP instance using `static_graph` (CSR); `source()`, `target()`, `is_tour()`
- `rcspp::sep::SeparationOracle` — Solver-independent cut separation (support graph, GH tree, parallel separators)
- `rcspp::preprocess::BoundPropagator` — Solver-independent domain propagation
- `rcspp::Model` — User-facing solver interface (HiGHS); `set_source()`/`set_target()` for paths, `set_depot()` for tours

## Coding Conventions

- C++23, GCC 14, `rcspp::` namespace (sub: `rcspp::sep::`, `rcspp::heuristic::`, `rcspp::preprocess::`)
- Tour vs s-t path: `source == target` → closed tour (degree 2); `source != target` → open path (degree 1 at endpoints)

## Dependencies

GCC 14, C++23, CMake 3.25+, TBB (`apt install libtbb-dev`), HiGHS + Catch2 (FetchContent)

## Workflow: Plan → Grind

Every task has two phases. Do not skip planning.

### 1. Plan (default)

When given a task, **start by planning**:

1. Investigate — read relevant code, docs, tests, issues
2. Propose an approach — what to change, where, and why
3. Discuss with the user — refine until aligned
4. Wait for explicit approval to proceed (e.g. "grind", "go", "do it")

### 2. Grind (on approval)

When the user says **"grind"** (or similar), execute autonomously:

1. Implement the change
2. Build — if it fails, read errors, fix, rebuild
3. Test — if tests fail, read failures, fix, retest
4. Repeat 1–3 until build and tests pass clean
5. Self-review: correctness, edge cases, performance
6. If review finds issues, go back to 1
7. Fullgate (see below) — branch, PR, sync, docs, push

Progress lives in files and git — not in your context window.

### When to Stop and Ask

Only pause the grind and ask a human when:
- A fix requires changing the public API or architecture
- You discover a bug in unrelated code you shouldn't touch
- You're stuck after multiple failed attempts at the same issue

Otherwise: keep going until build and tests pass.

### Fullgate

The final stage of a grind. Also runs standalone when user says **"fullgate"**:

1. **Branch** — create feature branch if needed
2. **PR** — create draft PR if none exists
3. **Sync** — pull latest main, merge into feature branch, resolve conflicts
4. **Tests** — add or update tests as needed
5. **Docs** — update README.md and docs as needed
6. **Push** — push branch, update PR description
7. **Review** — self-review: correctness, style, tests, performance
8. **Build** — `cmake --build build -j$(nproc)`
9. **Test** — `./build/rcspp_algo_tests && ./build/rcspp_tests`
10. **Push** — push any review fixes
11. **Finalize** — squash-merge PR, delete branch, pull main

### Claiming Work

- Add label `agent-wip` when you open or start working on an issue or PR
- Check for `agent-wip` before picking up work — never work on labeled items
- Remove `agent-wip` and close/merge when done

### Teams

When a task has independent sub-tasks, launch a team.
Each teammate runs in its own worktree (isolated repo copy).
Lead agent integrates: merge branches, resolve conflicts, run full build/test.
Do NOT use teams for sequential work.
