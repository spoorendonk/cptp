# rcspp-bac

`rcspp-bac` is a branch-and-cut solver for the Resource Constrained Shortest Path Problem (RCSPP) with profits, capacity, and optional customer visits.
Given a graph, edge costs, node profits, node demands, and capacity \(Q\), it finds a minimum-cost-minus-profit route as either a closed tour or an open s-t path.

The tour variant follows the CPTP formulation from Jepsen et al. (2014): *A branch-and-cut algorithm for the capacitated profitable tour problem*, Discrete Optimization 14, 78-96, https://doi.org/10.1016/j.disopt.2014.07.002.

## Mathematical Model

For graph \(G=(V,E)\), costs \(c_e\), profits \(p_i\), demands \(d_i\), and capacity \(Q\):

\[
\begin{aligned}
\min \quad & \sum_{e \in E} c_e x_e - \sum_{i \in V} p_i y_i \\\\
\text{s.t.} \quad
& \sum_{e \in \delta(i)} x_e =
\begin{cases}
y_i, & i \in \{s,t\} \\\\
2y_i, & i \in V \setminus \{s,t\}
\end{cases}
&& \forall i \in V \\\\
& \sum_{i \in V} d_i y_i \le Q \\\\
& y_s = 1,\; y_t = 1 \\\\
& x_e \in \{0,1\} && \forall e \in E \\\\
& y_i \in \{0,1\} && \forall i \in V
\end{aligned}
\]

Connectivity and capacity-defining inequalities are separated dynamically during branch-and-cut.

## Capabilities

- Delivers competitive performance on benchmark instance sets (results: [Benchmarks](docs/benchmarks.md)).
- Solves RCSPP instances as either closed tours or open s-t paths through one CLI/API workflow.
- Targets strong primal solutions and tighter dual bounds on difficult instances (how: [Primal Heuristic](docs/primal-heuristic.md), [Domain Propagator](docs/domain-propagator.md)).
- Strengthens LP relaxations during branch-and-cut with multiple cut families (how: [Cut Separation](docs/separation.md)).
- Reduces search space with preprocessing and route-bounding techniques (how: [Preprocessing](docs/preprocessing.md), [ng/DSSR Bounds](docs/ng-dssr-labelling.md)).
- Runs on HiGHS with a documented end-to-end solve pipeline (how: [Algorithms Overview](docs/algorithms.md)).

## Getting Started (CLI)

Build:

```bash
apt install libtbb-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Run:

```bash
./build/rcspp-solve <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]
```

Examples:

```bash
./build/rcspp-solve tests/data/tiny4.txt --time_limit 30 --output_flag false
./build/rcspp-solve tests/data/tiny4_path.txt
./build/rcspp-solve tests/data/tiny4.txt --source 0 --target 3
./build/rcspp_tests
pytest tests/python/test_solver.py
RCSPP_RUN_LOCAL_INSTALL_TEST=1 pytest tests/python/test_local_install.py
```

See:
- [CLI Usage + Parameters](docs/cli.md) (includes parameter table placeholder)
- [Instance Formats](docs/instance-formats.md)
- [C++ and Python API](docs/api.md)
