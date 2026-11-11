<!-- wcp: CLAUDE.md base=a4298469922d type=cpp-python -->
<!-- wcp-declined: 996da8f2cf7d Communication Style -->
<!-- wcp-declined: 0c56c7b46579 Code Navigation -->
<!-- wcp-declined: 058a3f5e0313 Issue Tracking -->
<!-- wcp-declined: a6e5ec6d33d8 CLAUDE.md Discipline -->
<!-- wcp-declined: 6e3ad83d1cb5 Complexity -->
<!-- wcp-declined: 448be1202c3e CMake -->
<!-- wcp-declined: e9feadcddcb5 Testing (GoogleTest) -->
<!-- wcp-declined: cc3df5ee540d Testing (pytest) -->
<!-- wcp-declined: 004504610dd6 Dependencies -->
<!-- wcp-declined: c4fc130d1235 Ownership and Lifetime -->
<!-- wcp-declined: cf8d0fc1773f Type Conversions -->
<!-- wcp-declined: 7e7de94f1273 Testing Bindings -->
<!-- wcp-declined: 9070a09f6181 LSP -->
# Project: cptp

Capacitated Profitable Tour Problem solver. Branch-and-cut solver for the CPTP following Jepsen et al. (2014). Also solves open s–t path variants. MIP backend is HiGHS with custom callbacks for user cut separation, domain propagation, and hyperplane branching.

## Build & Test

```clean
rm -rf build
```

```build
cmake -B build $(command -v ninja >/dev/null && echo -G Ninja) && cmake --build build -j$(nproc)
```

```test
ctest --test-dir build --output-on-failure -j$(nproc)
```

`pre-push` runs these three blocks as the local CI (`## Tooling and Gates`),
so they must need nothing but the build's own dependencies. The Python suite
runs inside ctest as `python_tests` when configure finds a Python 3.10+ with
pytest, numpy and its development headers (`.venv/bin` first, so `uv sync`
is enough): the bindings are then built too, and the tests import the
freshly built extension, not the venv's editable snapshot. Without one,
configure warns and the gate skips them. The lint tests (label `lint`) are
also ctest tests; `ctest --test-dir build -LE lint` skips them while
iterating.

Run a single C++ test by name:
```bash
./build/cptp_tests "test name pattern"
```

Python: each clone makes its own venv with `uv sync` (it compiles the
extension, HiGHS included, so on a shared machine run it under the build
budget, e.g. `CMAKE_BUILD_PARALLEL_LEVEL=4`). Never symlink another checkout's
`.venv`. pytest, ruff and mypy are in the `dev` dependency group, which
`uv sync` installs. Any `uv run` re-syncs, and so recompiles the extension,
after a change to `pyproject.toml`; the pre-commit hook runs `uv run` too, so
run `uv sync` under the build budget before committing such a change.
```bash
uv run pytest tests/python
uv run pytest tests/python/test_solver.py::test_name -x
```
The editable install snapshots the compiled extension, and uv decides whether
to rebuild from `pyproject.toml` alone, so after a change to C++ or
`python/bindings.cpp` run `uv sync --reinstall-package cptp` before the Python
suite.

Format C++ (every tracked source, as `pre-commit` checks them):
```bash
git ls-files -z '*.h' '*.cpp' | xargs -0 clang-format -i
```

The build patches HiGHS with custom user-cut and propagator callbacks from `third_party/highs_patch/`.

## Development Workflow

1. **Plan**: for non-trivial work, enter plan mode and align on approach before implementing. Small fixes can go straight to code.
2. **Implement**: write code. Nothing formats on save; `pre-commit` checks the format (`## Tooling and Gates`).
3. **Test**: run tests locally before considering work done.
4. **Review**: use the built-in code review (if available) or manually review changes before pushing.

## Git Workflow

Trunk-based development with linear history on main. Commit directly to main and push when local gates (review, build, test) pass.

Feature branches are optional for larger changes:
- Always branch from main. Run `git checkout main && git pull` before creating a branch.
- Never create a branch from another feature branch.
- Keep branches short-lived. Merge to main quickly.
- Use rebase or squash merge to maintain linear history — no merge commits on main.

## Tooling and Gates

Git hooks live in `.githooks/` and are tracked; they came from
[wcp](https://github.com/spoorendonk/wcp) (`wcp init`), see the sections below.
`core.hooksPath` is per-checkout and cannot be committed, so **a fresh clone has
no gates until `wcp init` (or `git config core.hooksPath .githooks`) is run in
it**.

- **pre-commit** — `clang-format --dry-run -Werror` over every tracked C++ file, with a clang-format of the major version `.clang-format-version` pins (another one skips the check, saying so). It checks, it does not reformat: run the format command above and commit again. For the Python half it runs `uv run ruff` and `uv run mypy` over the tree (both from the `dev` group), and skips each, saying so, when the venv does not have it. Without a `.venv`, `uv run` first builds one, compiling the extension inside the hook: run `uv sync` once beforehand. No tests at commit.
- **commit-msg** — Conventional Commits, subject ≤72 chars. `fixup!`/`squash!`/`amend!` pass through so `git rebase --autosquash` works, and `type!:` marks a breaking change (both local edits to wcp's hook).
- **pre-push** — refuses merge commits and unsquashed `fixup!`/`squash!`/`amend!` commits on the trunk (the latter a local edit to wcp's hook) and a HEAD without a review receipt (`## Self-Review`), then runs the clean/build/test blocks above in a fresh worktree of the pushed commit. That build fetches and compiles HiGHS from scratch, so a push takes minutes. A push changing only `*.md` files (except the root CLAUDE.md, whose blocks are the local CI's commands), the root LICENSE or the root `.gitignore` skips the build.
- **Lint as ctest tests** (`cmake/Lint.cmake`, label `lint`) — `clang_format` checks every tracked C++ file, and `clang_tidy:<path>` tidies one translation unit with warnings as errors, one test per source the build compiles (library, CLI, every test program, and `python/bindings.cpp` when the bindings are built). `.clang-tidy`'s `HeaderFilterRegex` reports on this repository's headers only. Every unit is tidied as C++23; `src/model/model.cpp` and `src/model/highs_bridge.cpp` can be because `third_party/highs_patch/apply_patch.cmake` moves `~HighsMipWorker()` out of line (clang rejects HiGHS's inline one under C++23). The lint tests need clang-format and clang-tidy **22.1.8**, the group's one toolchain version, pinned in `.clang-format-version` and in the `dev` group (`uv sync` puts them in `.venv/bin`, which `cmake/Lint.cmake` and the hooks search first); a missing one, or one of another major version, is a configure warning that drops its tests (a configure error when the `CI` environment variable is set). The version changes only for every repository of the group at once, through the hub. pre-push runs the suite without the `clang_tidy:` tests, then only those the push reaches, through its cached runner; the build block uses Ninja when it is installed, which gives the runner the dependency log it needs for that (without Ninja it runs every `clang_tidy:` test, uncached). A `build/` configured with another generator makes the build block abort with "generator does not match": run the `clean` block first. A changed header is tidied at push only through a one-line stub unit, not through its includers, and there is no scheduled full lint: **after a header change, run the whole lint yourself** (`ctest --test-dir build -L lint -j4`).
- **No GitHub Actions** — the local CI (pre-push; on a shared machine, run under `xk lock`, below) is the gate, and landing never waits on a GitHub run. The repository has no workflows; a GitHub workflow, if one is ever added, runs at most weekly (`schedule` plus `workflow_dispatch`), never on push or pull request. Nothing re-checks main after a push, so a `--no-verify` push or one from an ungated clone goes unnoticed unless someone checks `refs/reviews/` by hand: run `wcp init` in every clone.

The old hooks also auto-formatted and re-staged (with shfmt and ruff too),
applied safe clang-tidy fixes and ran the suite at commit, and ran shellcheck
as a push warning; those are gone. ruff (its complexity rule included) and
mypy run at commit.

`.clang-format` asks for `Standard: Latest`: the earlier `Standard: c++23` is
not a value clang-format accepts, so it refused to read the file and the tree
was never formatted to it. The whole tree was reformatted once, in the commit
listed in `.git-blame-ignore-revs`; run `git config blame.ignoreRevsFile
.git-blame-ignore-revs` to have `git blame` skip it.

On a shared machine, pin the local CI to the job budget: `git config
hooks.localCiCpus 0-3` makes the blocks' `-j$(nproc)` see four CPUs, and run the push itself under the machine lock (`xk lock -- git push …`).

## Hooks Are Files In This Repository

The git hooks live in `.githooks/` as ordinary shell scripts, committed like any other file. `wcp init` copies them there and points `core.hooksPath` at them; wcp installs no Claude Code hooks at all. The gate is those files and nothing else, so a clone without wcp installed is still gated — including the review receipt, which is a line of plain git the Self-Review section documents rather than anything wcp runs.

- **They are yours to read, and yours to edit.** Each hook carries a `# wcp:` line recording the digest of the bytes wcp wrote. That is what lets `wcp doctor` tell an edit you made from a version wcp has moved past: an untouched hook is refreshed silently, an edited one is left alone, and one that moved on both sides is reported rather than overwritten. Edit them when the project needs something different — and say so upstream if the default was wrong for everyone.
- **A clone is not gated until someone runs `wcp init` in it.** git clones neither a hooks directory nor local config, so the files arriving with the clone is not the same as git running them. `wcp doctor` reports an unset `core.hooksPath` as an error for exactly this reason.
- **They call no wcp.** A hook that shelled out to a missing tool would exit 127, which git reads as a refusal — so a machine without wcp could not commit, and could not commit the fix either. What the hooks do call is the project's own toolchain, and there they fail open and say so on stderr.
- **Never use `git push --no-verify` or `git commit --no-verify`** unless explicitly asked. A failing hook is a signal — fix the root cause, or change the hook deliberately and commit that.

## This Document Carries A Provenance Line

The comment at the top of this file is wcp's record of the bytes it wrote, in the same scheme the hooks use. `wcp doctor` reads it to tell an edit you made from a version wcp has moved past, and it is where this project's type is recorded — there is no `.wcp/` directory and no config file beside it.

**The prose is yours.** wcp writes this document once and never rewrites it unasked; edit it freely. When wcp adds or changes a section you have not edited, `wcp doctor` names it and `wcp init --interactive` offers it, one section at a time; answering `n` declines that version for good. Leave the top line in place, and the tooling stays able to tell your changes from its own.

## Running wcp

wcp is a command-line tool installed once per machine, outside any repository. The hooks never call it, so a machine without it is still gated.

```sh
uv tool install git+https://github.com/spoorendonk/wcp   # once per machine; @vX.Y.Z for a release
uv tool upgrade wcp                                       # to update an install without a tag
uv tool install --force git+https://github.com/spoorendonk/wcp@vX.Y.Z   # to move a tagged one
wcp --help                                                # the commands; `wcp <command> --help` for each one
```

- **Run `wcp doctor`** after pulling, after updating wcp, and before editing the hooks or this file. It changes nothing, and every finding names the command that fixes it. It exits 1 only when the gate does not work, e.g. `core.hooksPath` is unset, which is how every fresh clone starts until someone runs `wcp init` in it.
- **Files wcp wrote carry a `wcp:` provenance line**: `# wcp:` in each hook, and the comment at the top of this file. Edit the files when the project needs to, but leave that line in place. It is how `wcp doctor` tells your edits apart from a newer wcp version.
- **If `wcp` is not on PATH**, say so rather than installing it unasked. Installing it is a machine-wide change.

## Working Through Issues

A session's goal is often a list of issues (`wcp goal`). Work through them in the order given, one at a time:

1. **Issue check.** Read the issue and what it links. If it is clear and implementable, say so in one line and go on; this should be fast. If not, write a plan as a comment on the issue (plans go on issues, never in git) and work from it. If a later issue's plan changes an earlier one, re-plan then.
2. **Implement.**
3. **Gate and commit.** Run the project's gate (`## Build & Test`) and commit: a review reads commits, never uncommitted changes (`## Self-Review`).
4. **Review by one cold reviewer**: a fresh teammate or subagent running `/code-review` on the commits not yet pushed, never several reviewers in parallel. On findings, fix, gate, commit, and review again with a fresh reviewer. **At most two rounds**; the fix after the second is confirmed, not reviewed again (`## Self-Review`). Anything real left after that, or found out of scope, becomes its own issue. The bar for filing is a concrete failure scenario or a missing acceptance criterion; nits are fixed within the round or dropped.
5. **Record, push, close.** Record the review on a clean tree, push it with the commit (`## Self-Review`), and close the issue.

Without agent teams, implement in the session itself and review with a fresh subagent; the caps and the bar are the same.

**A hygiene goal** (`/goal Hygiene round …`, set by `wcp hygiene`) is a task, not a list of issues: work it from its note. Its findings are edits to existing issues or new issues labelled `tech-debt`; nothing is fixed inside the round, one cold reviewer checks the findings (at most two rounds), and the task ends with the summary comment its note specifies, whose first line is the marker wcp looks for.

## Self-Review

**Code you produced gets reviewed before it lands on main.** Run `/code-review`
on your own changes and act on what it finds; the review is not a formality you
report having done. The pre-push hook is the backstop, not the review.

**Then record the review, and push it with the commit**, or `pre-push` will
refuse the push:

```sh
git update-ref "refs/reviews/$(git rev-parse HEAD)" HEAD
git push --atomic origin HEAD "refs/reviews/$(git rev-parse HEAD)"
```

The receipt is a ref named by the reviewed commit, `refs/reviews/<sha>`, and it
travels with the push, so anyone can check by hand every commit that reaches the
trunk — including one pushed with `--no-verify`, or from a clone whose hooks
never ran. This repository runs no CI step for it (`## Tooling and Gates`).
`--atomic` lands the branch and its receipt together or not at all.

Plain git, on purpose. No wcp, no harness, nothing installed — so a reviewer on
any machine can write it, and a clone that has never heard of wcp is still able
to record a review and push. Every worktree of a clone sees the same receipts,
and each names one commit, so a review counts for that commit wherever it is
checked out and for nothing else.

**What it means, and what it does not.** It records that a review happened at
this commit — not what was reviewed. It carries exactly the trust of whoever ran
it, and it is a line anybody can write by hand; what it buys is that nobody can
push work no one has looked at *without noticing they are doing it*.

A later commit has no receipt, and the next push is refused again. That is
deliberate: a new commit is work nobody reviewed, including the commit that
fixes what the review found.

**The fix after the last round.** Reviews stop after two rounds, so the commit
that fixes the second round's findings gets no third review. That round's
reviewer reads the fix commit alone (`git show <sha>`) against its own
findings, and the receipt is recorded for the fix once it confirms they are
addressed. It looks for nothing new. What it does not confirm is taken out of
the fix, which makes a new commit: the same reviewer confirms that one before
its receipt is recorded, and what was taken out becomes its own issue.

**Record it only when the tree is clean.** With nothing new committed, a review
falls back to reading uncommitted changes — so a record written then names a
commit that does not contain what was read, and the push it clears is of work
nobody reviewed. Commit first, then review, then record. Untracked scratch files
do not matter; staged and unstaged changes to tracked files do.

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
- Type stubs: `python/cptp/_cptp.pyi` is generated by nanobind's stubgen and committed, so mypy (and a wheel's users, through `py.typed`) see the extension's real types. Every build of the bindings regenerates it in `build/python/cptp/`, and the `python_stub` ctest test fails when the committed copy differs. After changing `python/bindings.cpp`, refresh and commit it: `cmake --build build --target cptp_python_stub_update`. Don't edit it by hand; ruff skips it, mypy checks it.

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
