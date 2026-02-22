# Instance Formats

The solver auto-detects the format based on file extension and content.

## TSPLIB (`.vrp`, `.sppcc`)

Standard VRP / SPPRCLIB format. Detected by extension or by TSPLIB keywords in the file header.

```
NAME: example
TYPE: CVRP
DIMENSION: 5
CAPACITY: 10
EDGE_WEIGHT_TYPE: EXPLICIT
EDGE_WEIGHT_FORMAT: FULL_MATRIX
EDGE_WEIGHT_SECTION
 0  5  8 12  7
 5  0  6  9 11
 8  6  0  4 10
12  9  4  0  3
 7 11 10  3  0
DEMAND_SECTION
1 0
2 3
3 4
4 2
5 5
DEPOT_SECTION
1
-1
EOF
```

Optional `PROFIT_SECTION` is supported for instances with explicit node profits (used by the Roberti Set 3 instances). TSPLIB instances always use a tour formulation (depot from `DEPOT_SECTION`).

## Numeric (`.txt`)

Lightweight format for sparse graphs with optional s-t path support.

```
4 6               # num_nodes num_edges
0 1 10            # u v cost  (repeated for each edge)
0 2 8
0 3 12
1 2 6
1 3 7
2 3 5
7                 # capacity Q
0 3               # source target (optional; omit or use "0 0" for tour)
```

Node profits and demands follow the edge list when present. See `src/core/io.cpp` for the full specification.

When the `source target` line is omitted or `source == target`, the solver uses a tour formulation. When `source != target`, it uses an open s-t path formulation.
