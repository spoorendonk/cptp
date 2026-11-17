# cptp — Capacitated Profitable Tour Problem Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)
[![arXiv](https://img.shields.io/badge/arXiv-2607.04497-b31b1b.svg)](https://arxiv.org/abs/2607.04497)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844234.svg)](https://doi.org/10.5281/zenodo.20844234)

Branch-and-cut solver for the Capacitated Profitable Tour Problem (CPTP) following [Jepsen et al. (2014)](https://doi.org/10.1016/j.disopt.2014.08.001). Also solves open s–t path variants.

The implementation and the component study accompanying: S. Spoorendonk, [An open, reproducible branch-and-cut for the capacitated profitable tour problem: a component study](https://arxiv.org/abs/2607.04497), arXiv:2607.04497, 2026.

## Mathematical Model (Tour)

The tour variant uses `source = target = depot (0)`. The solver minimizes travel cost minus collected profits. `y_i` indicates whether node `i` is visited, and `x_e` indicates edge usage. Depot-incident edges are allowed to take value `2` to permit 2-node tours.

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

Here, `c_e` is edge cost, `p_i` node profit, `d^k_i` the demand of node `i` on capacity row `k`, `Q_k` that row's capacity, and `\delta(.)` denotes incident edges/cut boundary edges. An instance has m >= 1 capacity rows (the classical CPTP has m = 1), and a solution satisfies all of them. Connectivity (subtour elimination) inequalities are enforced via dynamic cut separation in branch-and-cut.

### Several capacity rows

Rows are numbered k = 1..m here; the code, its log and its option help call row 1 "row 0". A row is *bounded* when `Q_k < 1e17`; an unbounded row gives no capacity cuts or reachability fixing, and no hyperplanes unless it is row 1 (the labeling and the hyperplanes use row 1 even when it is unbounded). With m rows, each component uses by default:

| Component | Rows it uses by default | Option (default first) | Issue |
|---|---|---|---|
| Formulation | every row | — | [#90] |
| Reachability fixing (nodes that cannot fit) | every bounded row | — | [#90] |
| Primal heuristics | feasibility under every row; the LP-guided heuristic breaks ties by row 1's profit/demand | — | [#90] |
| Capacity cuts (RCI, multistar/GLM, RGLM) | every bounded row, each round | `capacity_cut_rows`: `all`, `most_violated`, `first` | [#91], measured on [#99] |
| Surrogate capacity cuts (RCI, RGLM on `sum_k d^k / Q_k`) | off | `capacity_cut_surrogate`: `false`, `true` | [#112], left off by [#113] |
| Labeling bounds: edge elimination, bounds propagation, reduced-cost fixing, all-pairs bounds (SPI) | row 1 | `labeling_rows`: `first`, `all`, `surrogate`, `vector` | [#101]–[#103], measured on [#104] |
| Edge elimination's join of the two labeling halves | by row 1's load | `labeling_join`: `load`, `min` | [#144] |
| 2-cycle elimination in the preprocessing labelings | off | `labeling_2cycle`: `false`, `true` | [#144] |
| Branching hyperplanes (`demand`, `clusters`) | row 1 and every bounded row k >= 2, each only when its demands are all integer; clusters are seeded by row 1's demands | `branch_hyper` (default `off`) | [#100] |

SEC and Comb cuts use no capacity row. Bounds from row 1 alone are valid, since a relaxation of one row is a relaxation of the intersection, but weaker than per-row ones: `labeling_rows=all` gains mainly at m = 10, at about 1.2x total time, and `vector` costs 3x at m = 10, with labelings of up to 972 s ([#104], [#115]). Cross-capacity inequalities, which combine rows, were studied in a literature write-up ([#105], on [#86]); surrogate capacity cuts closed none of the root gap at m = 2, 5 and 10 ([#113]). Further work on cross-row cuts is on [#146].

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

## Capabilities

- **MIP backend**: [HiGHS](https://highs.dev) with custom callbacks for user cut separation and domain propagation during branch-and-bound.
- **Cut separation**: SEC (Gomory-Hu), RCI, Multistar/GLM, RGLM, Comb, SPI (shortest-path inequalities — node-incompatibility cuts derived from all-pairs shortest-path bounds and a Held-Karp bound). SPI is a variant of the node-precedence / conflict-graph inequalities of García (2009), with roots in the shortest-path-bound preprocessing of Aneja et al. (1983); the Held-Karp bound is the only cptp-specific element.
- **Preprocessing / bounds**: Capacity-aware labeling with forward/backward bounds (on capacity row 1 by default; see [Several capacity rows](#several-capacity-rows)), optional all-pairs bounds, and bound-based edge/node elimination against the current upper bound.
- **Domain propagation**: sweep fixing (remove edges/nodes inconsistent with bound + UB checks), chain fixing (propagate implications from newly fixed edges), and Lagrangian reduced-cost fixing from LP reduced costs and bound checks.
- **Primal heuristics**: construction with 2-3 node seed routes, ILS neighborhoods (2-opt, relocate, swap, drop-add), and LP-guided in-tree ILS hooks on reduced graphs (threshold/RINS/neighborhood modes).
- **Branching**: hyperplane branching (`branch_hyper`, default `off`) adds one integer column z = a·y per hyperplane, linked by the row z - a·y = 0, and HiGHS branches on z as on any other column. The hyperplanes are fixed before the solve: pairs (Ryan-Foster-style customer pairs), clusters (small demand-seeded customer groups, weighed by each bounded capacity row's demands), demand (one per bounded capacity row) and cardinality (the number of selected customers); `all` adds every family. A capacity row whose demands are not all integer gives no demand or cluster hyperplane, since z would not be integer there.

## Build

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

## Run

```bash
./build/cptp-solve <instance> [--source <node>] [--target <node>] [--<option> <value> ...]

# Examples
./build/cptp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
./build/cptp-solve tests/data/tiny4.txt --source 0 --target 3
```

## Python

The Python package requires the C++ build dependencies (CMake, a C++23 compiler, libtbb-dev).

```bash
pip install .
```

```python
import cptp

# Load and solve an instance file
prob = cptp.load("instance.txt")
model = cptp.Model()
model.set_problem(prob)
result = model.solve([("time_limit", "60")])

# Or build from arrays
result = cptp.solve(
    num_nodes=4,
    edges=edges,
    edge_costs=costs,
    profits=profits,
    demands=demands,  # (n,), or (m, n) for m capacity rows
    capacity=10.0,  # a number, or m values
)

# Options of the multi-capacity components go to Model.solve
result = model.solve([("labeling_rows", "all"), ("capacity_cut_rows", "most_violated")])
```

`python -m cptp <instance> --<option> <value>` passes the same `--key value` options as `cptp-solve`, which `cptp-solve --help` lists; the Python CLI sets its own defaults for `time_limit` (600 s), `threads` (1) and `output_flag` (false, `--verbose` for true).

## Instance Formats

- **Numeric (.txt)**: native format — a first line `num_nodes num_arcs [m]`, nodes (id, profit, demand_1 .. demand_m), directed arcs (tail, head, cost; one arc per line), capacity_1 .. capacity_m on one line (when m = 1 it may be left out, for no limit, but only if the source/target line is left out too), optional source/target on one line, and nothing after it. m defaults to 1.
- **TSPLIB (.vrp)**: standard CVRP format with CAPACITY, DEMAND_SECTION, etc. With m rows, `CAPACITY : Q_1 .. Q_m`, and each DEMAND_SECTION line is `id d_1 .. d_m`.
- **SPPCC (.sppcc)**: SPPRCLIB column-generation subproblem format, with m rows as in TSPLIB.

`benchmarks/instances/multicap/` holds m = 2, 5 and 10 versions of the benchmark instances (`benchmarks/gen_multicap.py`). In Python, pass `demands` as an (m, n) array and `capacity` as m values.

## Benchmarks

Proven optimal within a 3600 s limit per instance (HiGHS v1.15.0), as reported in the paper ([arXiv:2607.04497](https://arxiv.org/abs/2607.04497)):

| Set | Instances | Solved | Rate |
|---|--:|--:|--:|
| SPPRCLIB (45–262 nodes) | 45 | 45 | 100% |
| Roberti Set 3 (45–200 nodes) | 31 | 24 | 77% |

Per-instance results are in `benchmarks/cptp.csv`; see `benchmarks/` for reproduction scripts.

**Per-component ablation.** `benchmarks/run_ablation.sh` re-runs the solver over the benchmark sets under several cut/fixing configurations (SEC-only, capacity cuts, Comb/RGLM, reduced-cost fixing on/off, SPI on/off) and writes one row per (config, instance) to `benchmarks/ablation.csv`.

### Known issues in earlier results

Two bugs, both present since the first commit (5dcf2b2, 2026-03-10), could
cut off an optimal solution and so report a suboptimal objective as optimal:

- **Reduced-cost fixing ([#114](https://github.com/spoorendonk/cptp/issues/114)).**
  The fixing test built an edge's bound on z_LP where the Lagrangian bound needs
  the base z_LP − rc·x\*, so it could fix out edges the optimum
  uses. RC fixing is on by default (`rc_fixing=adaptive`). Fixed in cf191b2.
- **Labeling bounds ([#111](https://github.com/spoorendonk/cptp/issues/111)).**
  The scalar labeling queued labels by their index in a node's list, and
  dominance pruning shifted those indices, so a live label could go
  unextended. The bounds then came out too high, which the edge/node
  elimination and the in-tree RC fixing trust. Fixed in 4cc1ee4.

Every result produced before cf191b2 (2026-11-12) may therefore hold wrong
optima. That includes the table above, every version of
`benchmarks/cptp.csv` and `benchmarks/ablation.csv`, the v0.1.0 and v0.2.0
releases, and the tables of the arXiv v1 / EJCO submission (2026-07-05).
They have not been rerun. Regenerate any result you rely on from cf191b2 or
later, the first commit with both fixes.

## Tests

Tests use multiple threads by default (parallel separation is exercised).

```bash
# C++ tests (Catch2-discovered; runs cptp_tests + cptp_tests_extra)
ctest --test-dir build --output-on-failure

# Optional direct test binaries
./build/cptp_tests
./build/cptp_tests_extra

# Python tests
pytest tests/python/test_solver.py
```

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

## License

MIT License — Copyright (c) 2026 Simon Spoorendonk

See [LICENSE](LICENSE) for the full text.
