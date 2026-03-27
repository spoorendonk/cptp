# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project: cptp — Capacitated Profitable Tour Problem Solver

Branch-and-cut solver for the CPTP following Jepsen et al. (2014). Also solves open s–t path variants. MIP backend is HiGHS with custom callbacks for user cut separation, domain propagation, and hyperplane branching.

## Build & Test

```build
cmake -B build && cmake --build build -j$(nproc)
```

```test
ctest --test-dir build --output-on-failure -j$(nproc) && pytest --tb=short -q
```

Run a single C++ test by name:
```bash
./build/cptp_tests "test name pattern"
```

Run a single Python test:
```bash
pytest tests/python/test_solver.py::test_name -x
```

Python editable install (for binding development):
```bash
pip install -e .
```

Format C++:
```bash
clang-format -i src/**/*.{h,cpp}
```

The build patches HiGHS with custom user-cut and propagator callbacks from `third_party/highs_patch/`.

## Development Workflow

1. **Plan**: for non-trivial work, enter plan mode and align on approach before implementing. Small fixes can go straight to code.
2. **Implement**: write code. Hooks auto-format, type-check, and flag complexity on save.
3. **Test**: run tests locally before considering work done.
4. **Review**: use the built-in code review (if available) or manually review changes before pushing.

## Git Workflow

Trunk-based development with linear history on main. Commit directly to main and push when local gates (review, build, test) pass.

Feature branches are optional for larger changes:
- Always branch from main. Run `git checkout main && git pull` before creating a branch.
- Never create a branch from another feature branch.
- Keep branches short-lived. Merge to main quickly.
- Use rebase or squash merge to maintain linear history — no merge commits on main.

## Commit Messages

Use Conventional Commits.

- Format: `type: description` or `type(scope): description`
- Types: `feat`, `fix`, `refactor`, `test`, `docs`, `style`, `perf`, `chore`, `build`, `ci`
- Subject line max 72 characters. Focus on **why**, not what.

## Architecture

The solver has a layered design:

- **`core/`** — Problem definition (`Problem`, `Edge`, `Solution`), I/O (multiple instance formats), and graph primitives (Gomory-Hu trees, Dinitz max-flow, static adjacency graph).
- **`sep/`** — Cut separators (SEC, RCI, Multistar/GLM, RGLM, Comb, SPI). Each implements the `Separator` interface. `SeparationOracle` orchestrates all separators per callback round.
- **`preprocess/`** — Bound propagation (capacity-aware labeling with forward/backward bounds), edge/node elimination against the current upper bound.
- **`heuristic/`** — LP-guided primal heuristic (construction + ILS neighborhoods: 2-opt, relocate, swap, drop-add).
- **`model/`** — `Model` is the user-facing API. `HiGHSBridge` wires the CPTP formulation, separators, propagator, and heuristic callback into HiGHS. This is the central integration point.
- **`cli/`** — Command-line entry point (`cptp-solve`).
- **`python/`** — nanobind bindings exposing `Model`, `Problem`, and `solve()` to Python.

Key data flow: `Problem` → `Model::solve()` → `HiGHSBridge` builds the IP formulation, installs separator/propagator/heuristic callbacks → HiGHS branch-and-cut → `SolveResult`.

## C++

- Target C++23. Use modern features (`std::expected`, concepts, ranges, `constexpr`).
- Style: Google-based, enforced by `.clang-format` and `.clang-tidy`.
- Use `#pragma once` for include guards.
- Minimize includes in headers. Forward-declare where possible.
- Use FetchContent for CMake dependencies. One `CMakeLists.txt` per directory with source files.
- Test files: `<module>_test.cpp` in `tests/`. Name tests descriptively: `TEST_F(SolverTest, ReturnsOptimalForFeasibleInput)`.

## Python

- Style: enforced by `ruff` (format + lint) and `mypy --strict`, configured in `pyproject.toml`.
- All functions must have full type annotations (mypy strict mode).
- Use built-in generics (`list[int]`, `dict[str, Any]`) and `|` union syntax.
- Test files: `test_<module>.py` in `tests/`. Use `conftest.py` for shared fixtures, `pytest.mark.parametrize` for data-driven tests.
- Pin dependencies with `>=` lower bounds in `pyproject.toml`.

## nanobind Bindings

- Bindings live in `python/`, separate from core C++ logic.
- C++ `camelCase` methods → Python `snake_case` via nanobind. Use `nb::arg("name")` for Python-friendly parameter names.
- Default: nanobind manages ownership. Use `nb::rv_policy::reference` only when C++ retains ownership and guarantees the object outlives Python references.
- Never return raw pointers without explicit lifetime annotation.
- Prefer returning by value or `std::shared_ptr`. Document ownership on each binding that transfers or shares it.
- Use automatic conversions for standard types (`std::string` ↔ `str`, `std::vector` ↔ `list`). Use `nb::ndarray` for NumPy interop.
- Test bindings from Python using pytest, not from C++.

## Reference Correctness

When implementing from papers, pseudocode, or open-source references:
- Match the reference algorithm exactly. Do not add early exits, iteration limits, size caps, or "optimization" shortcuts that change the algorithm's behavior.
- Only introduce heuristic approximations when explicitly asked.
- If the reference has edge cases or special handling, implement them — don't simplify them away.

## Plan Adherence

- Follow the agreed plan exactly. Do not take shortcuts, skip steps, or defer work without explicit permission.
- If you think the plan should change, stop and discuss — don't silently diverge.
- Implement everything specified in the plan. Do not leave TODO placeholders or stub implementations unless explicitly asked.

## Common Mistakes

### Behavioral
- **Don't add unrequested features.** If the user asked for X, deliver X — not X plus "helpful" extras.
- **Don't refactor surrounding code.** A bug fix doesn't need adjacent code cleaned up.
- **Don't skip tests to move faster.** Run the test suite even when changes seem trivial.
- **Don't silently change approach.** If something isn't working, say so — don't quietly try a different strategy.

### Technical
- **Don't invent APIs — verify they exist.** Check that functions, flags, and methods actually exist before using them.
- **Don't ignore type errors.** If mypy/clang-tidy flags something, fix the root cause — don't suppress or work around it.
- **Don't use deprecated patterns.** Check current docs, not training data. Prefer modern idioms.
- **Read before writing.** Always read a file before modifying it. Don't guess at existing code structure.
