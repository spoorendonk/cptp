# cptp — Capacitated Profitable Tour Problem Solver

Branch-and-cut solver for the Capacitated Profitable Tour Problem (CPTP) following [Jepsen et al. (2014)](https://doi.org/10.1016/j.disopt.2014.08.001). Also solves open s–t path variants.

## Capabilities

- **Cut separation**: SEC (Gomory-Hu), RCI, Multistar/GLM, RGLM, Comb, SPI (shortest-path inequalities — node-incompatibility cuts derived from all-pairs shortest-path bounds and Held-Karp DP; novel to this solver)
- **Primal heuristic**: deterministic warm-start incumbent finding (staged before HiGHS)
- **Domain propagation**: Lagrangian reduced-cost fixing (adaptive, periodic, root-only)
- **Preprocessing**: edge/node elimination, ng/DSSR route bounding
- **Concurrency**: deterministic staged preprocessing + deterministic callback scheduling
- **MIP backend**: [HiGHS](https://highs.dev) with hyperplane strong branching. Uses patched internal callbacks for user cut separation (cut pool injection via `HighsUserSeparator`) and domain propagation (Lagrangian fixing during B&B via propagator callback)

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

```bash
pip install cptp
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
| Roberti Set 3 (31 instances, 76–200 nodes, 300s) | 31 | 22 | 71% |

Single-thread, AMD Ryzen 9 3950X. See `benchmarks/` for full results and reproduction scripts.

## Tests

```bash
# C++ tests (Catch2-discovered; runs cptp_tests + cptp_tests_extra)
ctest --test-dir build --output-on-failure

# Optional direct test binaries
./build/cptp_tests
./build/cptp_tests_extra

# Python tests
pytest tests/python/test_solver.py
```

## License

MIT
