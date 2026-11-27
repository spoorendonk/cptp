# cptp — Capacitated Profitable Tour Problem Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)
[![arXiv](https://img.shields.io/badge/arXiv-2607.04497-b31b1b.svg)](https://arxiv.org/abs/2607.04497)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844234.svg)](https://doi.org/10.5281/zenodo.20844234)

cptp is a branch-and-cut solver for the Capacitated Profitable Tour Problem
(CPTP), following [Jepsen et al. (2014)](https://doi.org/10.1016/j.disopt.2014.08.001):
choose which customers to visit and in what order, so that travel cost minus
collected profit is smallest and the visited customers' demand fits the
vehicle's capacity. It also solves the open variant, a path from a source to a
target, and instances with several capacity rows, every one of which a solution
must satisfy. It is a C++23 library with a command-line solver and Python
bindings, on [HiGHS](https://highs.dev) as the MIP solver.

**Status:** complete, and not under active development. The accompanying
preprint, [arXiv:2607.04497](https://arxiv.org/abs/2607.04497), stays at v1.

What the solver brings to HiGHS's branch-and-bound:

- **Cuts**, separated in callbacks: subtour elimination (SEC, through
  Gomory–Hu trees), rounded capacity (RCI), multistar/GLM, rounded GLM (RGLM),
  comb, and shortest-path inequalities (SPI). SPI are node-incompatibility cuts
  from all-pairs shortest-path bounds and a Held–Karp bound; they are a variant
  of the node-precedence / conflict-graph inequalities of García (2009), with
  roots in the shortest-path-bound preprocessing of Aneja et al. (1983), and the
  Held–Karp bound is the only cptp-specific element.
- **Bounds and preprocessing**: a capacity-aware labeling gives forward and
  backward bounds, and optionally all-pairs bounds; edges and nodes that cannot
  be in a solution better than the current upper bound are eliminated.
- **Domain propagation**: sweep fixing (edges and nodes the bounds rule out),
  chain fixing (implications of newly fixed edges), and reduced-cost fixing
  from the LP's reduced costs.
- **Primal heuristics**: a construction from 2–3 node seed routes with iterated
  local search (2-opt, relocate, swap, drop-add), and an LP-guided search in the
  tree on reduced graphs (LP threshold, RINS, neighbourhood).
- **Hyperplane branching** (off by default): one integer column z = a·y per
  hyperplane, linked by the row z − a·y = 0, which HiGHS branches on as on any
  other column. The hyperplanes are fixed before the solve: pairs
  (Ryan–Foster-style customer pairs), clusters (small customer groups seeded by
  row 1's demands, weighed by the demands of row 1 and each bounded row),
  demand (one for row 1 and each bounded row) and cardinality (the number of
  selected customers); `all` adds every family. A capacity row whose demands
  are not all integer gives no demand or cluster hyperplane, since z would not
  be integer there.

## Quickstart

### Build

A C++23 compiler and CMake 3.25 or later; the build fetches HiGHS and patches
it, so the first build takes a few minutes.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
```

Set `--parallel` to the cores (and memory, about 2 GB a job) you can spare.

### Command line

```bash
./build/cptp-solve <instance> [--source <node>] [--target <node>] [--<option> <value> ...]
```

```console
$ ./build/cptp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
Instance: tiny4 (4 nodes, 6 edges, tour from depot 0)

Tour: 0 -> 1 -> 2 -> 0
Objective: -11  Bound: -11  Gap: 0%  Time: <seconds>s  Nodes: 1
...
$ ./build/cptp-solve tests/data/tiny4.txt --source 0 --target 3 --output_flag false
...
Path: 0 -> 1 -> 3
Objective: -13  Bound: -13  Gap: 0%  Time: <seconds>s  Nodes: 1
```

The objective is travel cost minus collected profit, so a negative one is a
gain. `--source` and `--target` turn a tour into an s–t path. The exit code is 0
when a solution was found.

### Python

`pip install` builds the C++ core with scikit-build-core and nanobind and
installs the `cptp` package and a `cptp` command, so it needs the same
compiler and CMake as the build above:

```bash
pip install .                                            # from a clone
pip install git+https://github.com/spoorendonk/cptp      # straight from GitHub
```

```python
import numpy as np
import cptp

# From an instance file
problem = cptp.load("tests/data/tiny4.txt")
model = cptp.Model()
model.set_problem(problem)
result = model.solve([("time_limit", "60"), ("output_flag", "false")])
print(result.status, result.objective, result.tour.tolist())

# From arrays: 4 nodes, depot 0, undirected edges
edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
costs = np.array([10.0, 8.0, 12.0, 6.0, 7.0, 5.0])
profits = np.array([0.0, 20.0, 15.0, 10.0])
demands = np.array([0.0, 3.0, 4.0, 2.0])
result = cptp.solve(4, edges, costs, profits, demands, capacity=7.0)
print(result.objective, result.tour.tolist())  # -11.0 [0, 1, 2, 0]

# Two capacity rows: demands is (m, n) and capacity has m values
demands = np.array([[0.0, 3.0, 4.0, 2.0], [0.0, 1.0, 3.0, 0.0]])
result = cptp.solve(4, edges, costs, profits, demands, capacity=[7.0, 3.0])
print(result.objective, result.tour.tolist())  # -1.0 [0, 1, 3, 0]
```

`cptp.solve` also takes `depot`, `source`, `target`, `time_limit`,
`num_threads` and `verbose`; any other option goes through `Model.solve`.
`cptp <instance> [--<option> <value> ...]` (or `python -m cptp ...`) is the
Python counterpart of `cptp-solve`.

## Instance formats

`cptp.load` and `cptp-solve` read `.vrp`, `.tsp` and `.sppcc` files as TSPLIB,
and any other file as TSPLIB when its content is, else as numeric:

- **Numeric (`.txt`)**, the native format: a first line `num_nodes num_arcs [m]`
  (m capacity rows, default 1); one line per node, `id profit demand_1 ..
  demand_m`; one line per arc, `tail head cost`; then the capacities
  `Q_1 .. Q_m` on one line; then, optionally, `source target` on one line, and
  nothing after it. The graph is undirected: an arc with tail < head gives the
  edge at its cost, and one with tail > head is read and dropped without a
  warning, so a file lists each edge either once with tail < head or both
  ways. With m = 1 the capacity line may be left
  out, for no limit, but only when the source/target line is left out too.
  [`tests/data/tiny4.txt`](tests/data/tiny4.txt) is a tour and
  [`tests/data/tiny4_path.txt`](tests/data/tiny4_path.txt) a path.
- **TSPLIB (`.vrp`)**: the CVRP format, with `CAPACITY` and `DEMAND_SECTION`.
  With m rows, `CAPACITY : Q_1 .. Q_m`, and each `DEMAND_SECTION` line is
  `id d_1 .. d_m`.
- **SPPCC (`.sppcc`)**: the SPPRCLIB column-generation subproblem format, with
  m rows as in TSPLIB.

[`benchmarks/instances/`](benchmarks/instances) holds the SPPRCLIB and Roberti
benchmark sets, and `multicap/` versions of them with m = 2, 5 and 10 capacity
rows (made by [`benchmarks/gen_multicap.py`](benchmarks/gen_multicap.py)).

## Options

Options are `--key value` pairs on `cptp-solve` and `python -m cptp`, and
`(key, value)` string pairs in `Model.solve`. **[docs/options.md](docs/options.md)
lists them with their defaults**: `cptp-solve --help` word for word, and the
defaults and tuning knobs the help leaves out. The ones most often set:

| Option | What it does |
|---|---|
| `time_limit` | Wall-clock limit in seconds, preprocessing included |
| `threads` | Threads for HiGHS and for cptp's own parallel work |
| `random_seed` | Seed for HiGHS and the heuristics |
| `output_flag` | Log to the console (`true` on `cptp-solve` and in `Model.solve`, `false` in `python -m cptp` and `cptp.solve`) |
| `cutoff` | A known upper bound, to prove it optimal |
| `parallel` | HiGHS's parallel tree search; `off` for a reproducible solve with more than one thread |

A key cptp does not know goes to HiGHS (`cptp-solve --highs_help` lists
HiGHS's options); presolve stays off. With more than one thread, node counts and,
among tied optima, the tour can differ between runs with the same seed;
`parallel off` or `threads 1` gives the same run every time.

## Model

The tour variant has `source = target = depot (0)`. `y_i` says whether node `i`
is visited and `x_e` whether edge `e` is used; edges at the depot may take the
value 2, which allows a tour to one customer and back.

$$
\begin{aligned}
\min \quad & \sum_{e \in E} c_e x_e - \sum_{i \in V} p_i y_i \\
\text{s.t.} \quad
& \sum_{e \in \delta(i)} x_e = 2 y_i && \forall i \in V \\
& y_0 = 1 \\
& \sum_{i \in V} d^k_i y_i \le Q_k && k = 1, \dots, m \\
& \sum_{e \in \delta(S)} x_e \ge 2 y_i && \forall S \subset V \setminus \{0\}, \ \forall i \in S \\
& x_e \in \{0,1\} && \forall e \in E \setminus \delta(0) \\
& x_e \in \{0,1,2\} && \forall e \in \delta(0) \\
& y_i \in \{0,1\} && \forall i \in V
\end{aligned}
$$

`c_e` is the edge cost, `p_i` the node profit, `d^k_i` the demand of node `i` on
capacity row `k`, `Q_k` that row's capacity, and `δ(S)` the edges with one end
in `S`. The classical CPTP has m = 1. The connectivity (subtour elimination)
constraints are separated as cuts during branch-and-cut, not stated up front.
The s–t path variant visits its source and target (`y = 1` there) with degree
1 instead of 2, and allows no value 2.

### Several capacity rows

Rows are numbered k = 1..m here; the code, its log and `--help` call row 1
"row 0". A row is *bounded* when `Q_k < 1e17`. An unbounded row gives no
capacity cuts and no reachability fixing, and no hyperplanes unless it is row 1
(the labeling and the hyperplanes use row 1 even when it is unbounded). A row
with `Q_k = 0` gives no capacity cuts either, nor a hyperplane for k ≥ 2;
reachability fixing removes its customers of positive demand.

With m rows, each component uses by default:

| Component | Rows it uses by default | Option (default first) |
|---|---|---|
| Formulation | every row | — |
| Reachability fixing (nodes that cannot fit) | every bounded row | — |
| Primal heuristics | feasibility under every row; the LP-guided heuristic breaks ties by row 1's profit/demand | — |
| Capacity cuts (RCI, multistar/GLM, RGLM) | every bounded row, each round | `capacity_cut_rows`: `all`, `most_violated`, `first` |
| Surrogate capacity cuts (RCI, RGLM on `sum_k d^k / Q_k`) | off | `capacity_cut_surrogate`: `false`, `true` |
| Multistar/GLM also on each target's largest depot–target min-cut side connected to the target, not only the Gomory–Hu tree's smallest (any m) | off | `multistar_largest_side`: `false`, `true` |
| Cross-row lifted cover w-GLMs: per bounded row a cover over the min-cut sides and V minus the depot, lifted against every bounded row | off | `enable_xcover`: `false`, `true` |
| Labeling bounds: edge elimination, bounds propagation, reduced-cost fixing, all-pairs bounds (SPI) | row 1 | `labeling_rows`: `first`, `all`, `surrogate`, `vector` |
| Edge elimination's join of the two labeling halves | by row 1's load | `labeling_join`: `load`, `min` |
| 2-cycle elimination in the preprocessing single-row labelings (not `labeling_rows=vector`'s vector labeling) | off | `labeling_2cycle`: `false`, `true` |
| Branching hyperplanes (`demand`, `clusters`) | row 1 and every bounded row k ≥ 2, each only when its demands are all integer; clusters are seeded by row 1's demands | `branch_hyper` (default `off`) |

SEC and comb cuts use no capacity row. Bounds from row 1 alone are valid, since
a relaxation of one row is a relaxation of the intersection, but weaker than per-row
ones: `labeling_rows=all` gains mainly at m = 10, at about 1.2x the total time,
and `vector` costs 3x at m = 10, with labelings of up to 972 s. Surrogate
capacity cuts closed none of the root gap at m = 2, 5 and 10, which is why they
are off. [docs/development.md](docs/development.md#where-the-defaults-come-from)
links the measurements behind each default.

## Results

### Multi-capacity smoke test

*Not yet run.* A smoke test of the solver against PathWyse 1.1.0 (its
`PWDefault` and `PWBucket` algorithms) on several capacity rows is scheduled on
[#147](https://github.com/spoorendonk/cptp/issues/147); its table goes here:
solved counts and the shifted geometric mean time per solver for m = 1, 2, 5
and 10. Read it as a smoke test, not a benchmark: 10 base instances (six
SPPRCLIB, four Roberti), a 600 s limit, and runs made concurrently on one
machine.

### One capacity row: the comparison with labeling

The preprint compares the branch-and-cut with PathWyse, a labeling
(dynamic-programming) solver, on SPPRCLIB and Roberti Set 3: PathWyse solved
71 of 76 instances and the branch-and-cut 69, with shifted geometric mean times
of 3.66 s and 31.10 s respectively; the branch-and-cut solved all five
instances PathWyse did not. That comparison used a PathWyse of about version 0.2.0, against the
branch-and-cut before the fixes below. PathWyse 1.1.0 reports itself 26.7x
faster on SPPRCLIB (`PWBucket` 1.1.0 against `PWDefault` 1.0.0), so on one
capacity row, labeling is the faster method.
[#147](https://github.com/spoorendonk/cptp/issues/147) pins PathWyse 1.1.0 and
holds the smoke test above.

### Known issues in earlier results

Two bugs, both present since the first commit (5dcf2b2, 2026-03-10), could cut
off an optimal solution and so report a suboptimal objective as optimal:

- **Reduced-cost fixing ([#114](https://github.com/spoorendonk/cptp/issues/114)).**
  The fixing test built an edge's bound on z_LP where the Lagrangian bound needs
  the base z_LP − rc·x\*, so it could fix out edges the optimum uses. RC fixing
  is on by default (`rc_fixing=adaptive`). Fixed in cf191b2.
- **Labeling bounds ([#111](https://github.com/spoorendonk/cptp/issues/111)).**
  The scalar labeling queued labels by their index in a node's list, and
  dominance pruning shifted those indices, so a live label could go unextended.
  The bounds then came out too high, which the edge/node elimination and the
  in-tree RC fixing trust. Fixed in 4cc1ee4.

Every result produced before cf191b2 (2026-11-12) may therefore hold wrong
optima. That includes the table below, every version of
[`benchmarks/cptp.csv`](benchmarks/cptp.csv) and
[`benchmarks/ablation.csv`](benchmarks/ablation.csv), the v0.1.0 and v0.2.0
releases, and the tables of the preprint (arXiv v1, 2026-07-05). They have not
been rerun. Regenerate any result you rely on from cf191b2 or later, the first
commit with both fixes.

A third bug left SPI without its all-pairs bounds
([#156](https://github.com/spoorendonk/cptp/issues/156)): the propagator moved
them out of the bridge before the solve, so the SPI separator always saw an
empty matrix and returned no cut. **Every result produced with
`enable_spi=true` before the fix ran without SPI**, the `spi_on` rows of
`benchmarks/ablation.csv` included (their `spi_cuts` are 0, or empty where the
run left no stats); `spi_on` vs `spi_off` measured nothing about SPI.

**Earlier results, before the fixes (may hold wrong optima).** The preprint
reports these instances proven optimal within 3600 s each, with HiGHS v1.15.0:

| Set | Instances | Solved | Rate |
|---|--:|--:|--:|
| SPPRCLIB (45–262 nodes) | 45 | 45 | 100% |
| Roberti Set 3 (45–200 nodes) | 31 | 24 | 77% |

The per-instance rows are in `benchmarks/cptp.csv`, and the per-component
ablation (SEC only, capacity cuts, comb/RGLM, reduced-cost fixing and SPI on
and off) in `benchmarks/ablation.csv`. [`benchmarks/`](benchmarks) holds the
scripts that produce them: `run_benchmarks.sh`, `run_ablation.sh`, and
`run_pathwyse_reference.sh` with `compare_pathwyse.sh` for the PathWyse
comparison.

## Citing

Cite the paper for the result and the archived release for the code.
[`CITATION.cff`](CITATION.cff) carries both in machine-readable form — GitHub's
"Cite this repository" button reads it.

The paper is **arXiv:2607.04497**, <https://arxiv.org/abs/2607.04497> — a
component study of this solver: which cuts, bounds, propagation and branching
rules actually carry the branch-and-cut, measured by ablation.

```bibtex
@misc{spoorendonk2026cptp,
  title        = {An open, reproducible branch-and-cut for the capacitated
                  profitable tour problem: a component study},
  author       = {Spoorendonk, Simon},
  year         = {2026},
  eprint       = {2607.04497},
  archivePrefix= {arXiv},
  primaryClass = {math.OC},
  url          = {https://arxiv.org/abs/2607.04497}
}
```

The archived code is on Zenodo. **v0.2.0 is [10.5281/zenodo.21205088](https://doi.org/10.5281/zenodo.21205088)**; the
concept DOI [10.5281/zenodo.20844234](https://doi.org/10.5281/zenodo.20844234) always resolves to the latest version.
Cite the version DOI when reproducibility matters — it pins one tree and one set
of benchmark results. The earlier v0.1.0 is
[10.5281/zenodo.20844235](https://doi.org/10.5281/zenodo.20844235).

```bibtex
@software{spoorendonk2026cptpsoftware,
  title     = {cptp: a C++23 branch-and-cut solver for the Capacitated
               Profitable Tour Problem},
  author    = {Spoorendonk, Simon},
  year      = {2026},
  version   = {0.2.0},
  doi       = {10.5281/zenodo.21205088},
  url       = {https://doi.org/10.5281/zenodo.21205088}
}
```

Cite a tag or a DOI, never `main`: the benchmark tables move with the code.

## Development

Tests, the local gate, parallel-search details, ThreadSanitizer and the
measurements behind the defaults are in [docs/development.md](docs/development.md).
Release notes are on the [GitHub releases](https://github.com/spoorendonk/cptp/releases).

## License

MIT License — Copyright (c) 2026 Simon Spoorendonk

See [LICENSE](LICENSE) for the full text.
