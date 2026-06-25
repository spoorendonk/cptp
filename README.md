# cptp — Capacitated Profitable Tour Problem Solver

[![CI](https://github.com/spoorendonk/cptp/actions/workflows/ci.yml/badge.svg)](https://github.com/spoorendonk/cptp/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)
<!-- DOI badge added in the post-release backfill commit:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Branch-and-cut solver for the Capacitated Profitable Tour Problem (CPTP) following [Jepsen et al. (2014)](https://doi.org/10.1016/j.disopt.2014.08.001). Also solves open s–t path variants.

## Mathematical Model (Tour)

The tour variant uses `source = target = depot (0)`. The solver minimizes travel cost minus collected profits. `y_i` indicates whether node `i` is visited, and `x_e` indicates edge usage. Depot-incident edges are allowed to take value `2` to permit 2-node tours.

$$
\begin{aligned}
\min \quad & \sum_{e \in E} c_e x_e - \sum_{i \in V} p_i y_i \\
\text{s.t.} \quad
& \sum_{e \in \delta(i)} x_e = 2 y_i && \forall i \in V \\
& y_0 = 1 \\
& \sum_{i \in V} d_i y_i \le Q \\
& \sum_{e \in \delta(S)} x_e \ge 2 y_i && \forall S \subset V \setminus \{0\}, \ \forall i \in S \\
& x_e \in \{0,1\} && \forall e \in E \setminus \delta(0) \\
& x_e \in \{0,1,2\} && \forall e \in \delta(0) \\
& y_i \in \{0,1\} && \forall i \in V
\end{aligned}
$$

Here, `c_e` is edge cost, `p_i` node profit, `d_i` node demand, `Q` vehicle capacity, and `\delta(.)` denotes incident edges/cut boundary edges. Connectivity (subtour elimination) inequalities are enforced via dynamic cut separation in branch-and-cut.

## Capabilities

- **MIP backend**: [HiGHS](https://highs.dev) with custom callbacks for user cut separation, domain propagation, and hyperplane branching during branch-and-bound.
- **Cut separation**: SEC (Gomory-Hu), RCI, Multistar/GLM, RGLM, Comb, SPI (shortest-path inequalities — node-incompatibility cuts derived from all-pairs shortest-path bounds and a Held-Karp bound). SPI is a variant of the node-precedence / conflict-graph inequalities of García (2009), with roots in the shortest-path-bound preprocessing of Aneja et al. (1983); the Held-Karp bound is the only cptp-specific element.
- **Preprocessing / bounds**: Capacity-aware labeling with forward/backward bounds, optional all-pairs bounds, and bound-based edge/node elimination against the current upper bound.
- **Domain propagation**: sweep fixing (remove edges/nodes inconsistent with bound + UB checks), chain fixing (propagate implications from newly fixed edges), and Lagrangian reduced-cost fixing from LP reduced costs and bound checks.
- **Primal heuristics**: construction with 2-3 node seed routes, ILS neighborhoods (2-opt, relocate, swap, drop-add), and LP-guided in-tree ILS hooks on reduced graphs (threshold/RINS/neighborhood modes).
- **Branching**: hyperplane branching supports pairs (Ryan-Foster-style customer pairs), clusters (small demand-seeded customer groups), demand (weighted sum of selected demands), and cardinality (number of selected customers).

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
    num_nodes=4, edges=edges, edge_costs=costs,
    profits=profits, demands=demands, capacity=10.0,
)
```

## Instance Formats

- **Numeric (.txt)**: native format — nodes (id, profit, demand), directed arcs (tail, head, cost), capacity, optional source/target
- **TSPLIB (.vrp)**: standard CVRP format with CAPACITY, DEMAND_SECTION, etc.
- **SPPCC (.sppcc)**: SPPRCLIB column-generation subproblem format

## Benchmarks

| Set | Instances | Solved | Rate |
|---|--:|--:|--:|
| SPPRCLIB (45 instances, 45–262 nodes) | 45 | 45 | 100% |
| Roberti Set 3 (31 instances, 76–200 nodes, 3600s) | 31 | 27 | 87% |

See `benchmarks/` for full results and reproduction scripts.

**Per-component ablation.** `benchmarks/run_ablation.sh` re-runs the solver over the benchmark sets under several cut/fixing configurations (SEC-only, capacity cuts, Comb/RGLM, reduced-cost fixing on/off, SPI on/off) and writes one row per (config, instance) to `benchmarks/ablation.csv`.

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

## Citation

If you use this software, please cite the archived release via its Zenodo DOI.
Machine-readable metadata lives in [`CITATION.cff`](CITATION.cff), from which
GitHub renders a "Cite this repository" button.

<!-- Filled in by the post-release backfill commit, once Zenodo has minted the DOI:
     add the DOI badge to the header badge block above, then replace this comment
     with the citation, e.g.:

> Spoorendonk, S. (2026). *cptp* (v0.1.0). Zenodo. <https://doi.org/10.5281/zenodo.XXXXXXX>

The concept DOI `10.5281/zenodo.XXXXXXX` always resolves to the latest release;
the v0.1.0 version DOI is `10.5281/zenodo.YYYYYYY`. -->

## License

MIT License — Copyright (c) 2026 Simon Spoorendonk

See [LICENSE](LICENSE) for the full text.
