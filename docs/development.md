# Development

How to build, test and change cptp. The README is for using the solver; this
page is for working on it.

## Build and gate

The local gate is the three blocks under `## Build & Test` in
[CLAUDE.md](https://github.com/spoorendonk/cptp/blob/main/CLAUDE.md): clean, build (Ninja when installed), and the full
ctest suite. The `pre-push` hook in `.githooks/` runs them in a fresh worktree
of the pushed commit, after checking for a review receipt; a fresh clone has no
hooks until `wcp init` (or `git config core.hooksPath .githooks`) is run in it.
CLAUDE.md also covers formatting, the lint tests (label `lint`), the Python
venv (`uv sync`) and the review receipt.

## Tests

The tests use several threads by default: separation runs in parallel at the
root, before HiGHS's tree workers start, and in `parallel=off` runs.

```bash
# C++ tests (Catch2-discovered; runs cptp_tests and cptp_tests_extra)
ctest --test-dir build --output-on-failure

# The test binaries directly, with an optional name pattern
./build/cptp_tests
./build/cptp_tests_extra

# Python tests (also run inside ctest as python_tests when a venv is found)
uv run pytest tests/python
```

## Generated files

Two committed files are generated from the build, and a ctest test fails when
the committed copy differs:

| File | Test | Refresh with |
|---|---|---|
| `python/cptp/_cptp.pyi` (type stubs) | `python_stub` | `cmake --build build --target cptp_python_stub_update` |
| `docs/options.md` (the `cptp-solve --help` block) | `options_doc` | `cmake --build build --target cptp_options_doc_update` |

The options target rewrites only the block between the two marker comments in
`docs/options.md`; the prose around it is edited by hand.

## Documentation site

The site at <https://spoorendonk.dk/cptp/> is built with Sphinx from `docs/`:
`docs/index.md` pulls in the README, and next to it are this page, the
options, a C++ reference (Doxygen XML of the user-facing headers, rendered by
Breathe) and a Python reference (sphinx-autoapi, read statically from
`python/cptp/` and the committed `_cptp.pyi`, so the extension is not
compiled). It also writes `llms.txt` and `llms-full.txt`. Build it with
`doxygen` on PATH:

```bash
uv run --only-group docs sphinx-build -W -b html docs build/docs/html
```

`-W` turns warnings into errors; the build is warning-free. `--only-group docs`
installs the Sphinx packages into the venv without building cptp.

The gate keeps it warning-free: the ctest test `docs_build` (label `docs`)
runs the same build, into the build directory's `docs/html`, when configure finds `doxygen` and
the venv's Python has the `docs` group, and a push changing the README or a
page under `docs/` runs the gate. Without either, configure warns and the test
is not registered. `uv sync --group docs` keeps the group in the venv beside
`dev`; a plain `uv sync` removes it again. The `docs`
workflow (`.github/workflows/docs.yml`) builds and deploys the site to GitHub
Pages; it runs on a `v*` tag or by hand (`gh workflow run docs.yml`), never on
a push to a branch or a pull request.

## Parallel tree search

HiGHS's parallel tree search is on by default (HiGHS's own `parallel=choose`)
whenever a solve has more than one thread; `python -m cptp` defaults to
`--threads 1`, so it needs `--threads N` to use it. `parallel=off` runs one
branch-and-bound tree, with cptp's separators and labeling parallel inside it
instead. On a smoke set (six SPPRCLIB instances and their two-row multicap
versions, three seeds, `threads=4`, Mac mini) the parallel search took a median
0.81 of the serial tree's time (IQR 0.62–0.91) and reached the same optima
([#154]). The callbacks are safe under it ([#152], [#153]).

With more than one thread, its result depends on timing: node counts,
statistics and, among tied optima, the tour can differ between runs with the
same seed; `parallel=off` or `threads=1` gives a reproducible solve. Its
acceptance matrix (threads 2/4/8 × three seeds) is hidden from the default run:
`./build/cptp_tests "[parallel-matrix]"`.

## ThreadSanitizer

Configure a separate build with `-fsanitize=thread` in the global flags, so the
fetched HiGHS is instrumented too, and run the hidden stress test, every
instance under real parallel tree search at four threads (the
simulated-concurrency tests run on one thread, so TSan finds no race in them):

```bash
cmake -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_C_FLAGS=-fsanitize=thread -DCMAKE_CXX_FLAGS=-fsanitize=thread \
  -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=thread \
  -DCMAKE_SHARED_LINKER_FLAGS=-fsanitize=thread
cmake --build build-tsan --target cptp_tests
TSAN_OPTIONS=symbolize=0 ./build-tsan/cptp_tests "[parallel-stress]"
```

`symbolize=0` because on macOS TSan's symbolizer (`atos`) can hang on a report;
symbolize the addresses afterwards with `atos -o` if one shows up. No
HiGHS-internal suppressions are needed (2026-11-18, HiGHS v1.16).

## Where the defaults come from

The README's [Several capacity rows](https://github.com/spoorendonk/cptp#several-capacity-rows) table
gives each component's default; these issues hold the work and the
measurements behind them:

| Component | Issues |
|---|---|
| Formulation, reachability fixing, primal heuristics on every row | [#90] |
| Capacity cuts per row (`capacity_cut_rows`) | [#91], measured on [#99] |
| Surrogate capacity cuts (`capacity_cut_surrogate`) | [#112], left off by [#113] |
| Multistar/GLM on the largest min-cut side (`multistar_largest_side`) | [#168] |
| Cross-row lifted cover w-GLMs (`enable_xcover`) | [#166] |
| Labeling bounds per row (`labeling_rows`) | [#101]–[#103], measured on [#104] and [#115] |
| Labeling join and 2-cycle elimination (`labeling_join`, `labeling_2cycle`) | [#144] |
| Branching hyperplanes per row (`branch_hyper`) | [#100] |
| Root cut loop limits (`root_sepa_rounds`, `root_sepa_stall`) | [#150] |

Cross-capacity inequalities, which combine rows, were studied in a literature
write-up ([#105], on [#86]); surrogate capacity cuts closed none of the root gap
at m = 2, 5 and 10 ([#113]). Further work on cross-row cuts is on [#146].

[#86]: https://github.com/spoorendonk/cptp/issues/86
[#90]: https://github.com/spoorendonk/cptp/issues/90
[#91]: https://github.com/spoorendonk/cptp/issues/91
[#99]: https://github.com/spoorendonk/cptp/issues/99
[#100]: https://github.com/spoorendonk/cptp/issues/100
[#101]: https://github.com/spoorendonk/cptp/issues/101
[#103]: https://github.com/spoorendonk/cptp/issues/103
[#104]: https://github.com/spoorendonk/cptp/issues/104
[#105]: https://github.com/spoorendonk/cptp/issues/105
[#112]: https://github.com/spoorendonk/cptp/issues/112
[#113]: https://github.com/spoorendonk/cptp/issues/113
[#115]: https://github.com/spoorendonk/cptp/issues/115
[#144]: https://github.com/spoorendonk/cptp/issues/144
[#146]: https://github.com/spoorendonk/cptp/issues/146
[#150]: https://github.com/spoorendonk/cptp/issues/150
[#152]: https://github.com/spoorendonk/cptp/issues/152
[#153]: https://github.com/spoorendonk/cptp/issues/153
[#154]: https://github.com/spoorendonk/cptp/issues/154
[#166]: https://github.com/spoorendonk/cptp/issues/166
[#168]: https://github.com/spoorendonk/cptp/issues/168
