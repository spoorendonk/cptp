This is the version behind the results in the accompanying paper (arXiv / EJCO).

## Highlights

- **HiGHS 1.15.0.** Upgraded the MIP backend from 1.14.0 to 1.15.0 and ported the
  user-cut / propagator patch (user-cut separation, domain propagation, hyperplane
  branching, root-separation anchors) to the new release.
- **Correctness fix.** A warm-start path that could yield **wrong optima under HiGHS
  1.15** is fixed by routing the incumbent through the solution callback. Objectives
  now match the DP/labeling reference (PathWyse) on every commonly solved instance.
- **Benchmark instrumentation.** `cptp-solve` now reports per-phase timing
  (separation / propagation / heuristic / reduced-cost fixing), propagation and
  reduced-cost-fixing frequencies, and an LP-iteration breakdown
  (strong-branching / separation / heuristic); these are captured in `cptp.csv` /
  `ablation.csv`. Benchmark sweeps are chunked and resumable.
- **Results regenerated** on the fixed HiGHS 1.15 build.

## Changes since v0.1.0

Backend & correctness
- Upgrade to HiGHS 1.15.0 and port the user-cut/propagator patch.
- Restore propagator + root-separation anchors for HiGHS 1.15.
- Route warm-start via the callback to fix wrong optima on HiGHS 1.15.

Benchmarks & instrumentation
- Surface per-phase stats and chunked, resumable sweeps.
- Parse LP-iteration columns in `run_benchmarks.sh`.
- Regenerate `ablation.csv` with wall-clock phase data.
- Regenerate `cptp.csv` / `ablation.csv` on the fixed HiGHS 1.15 build.

Docs, tests & fixes
- Add Zenodo DOI to citation metadata.
- Correct Roberti node range and Jepsen citation year.
- Cover `setSolution` warm-start fallback and per-phase stats with tests.
- Address review nits (regression test, guard, comment).
