# Cut Separation

File: `src/sep/` directory

## Overview

Five families of cutting planes are separated dynamically. All separators are
**solver-independent** — they implement `rcspp::sep::Separator` and receive a
`SeparationContext` containing the LP relaxation solution, problem data, and a
shared Gomory-Hu tree. The `SeparationOracle` class bundles all separators into
a single call usable from any MIP solver's cut callback.

All cuts use formulations valid for RCSPP with optional vertices (y variables),
originally developed for the CPTP (Jepsen et al. 2014). Standard CVRP cuts are
**invalid** here because they assume all nodes are visited.

## Gomory-Hu Tree

File: `src/core/gomory_hu.h`

Implements **Gusfield's simplified algorithm** to compute all pairwise min-cuts
with only n-1 max-flow calls (instead of O(n^2)). The tree is built once per
separation round from the fractional support graph and shared across all separators.

Each max-flow uses the **Dinitz algorithm** (O(V^2 E)) on the support graph
(`src/core/dinitz.h`).

Two query methods:

- `min_cut(target)` — walks the tree path from target to root (depot), returns
  the bottleneck cut partition and value. Used by SEC, Multistar, RGLM.
- `all_cuts_on_path(target)` — returns all cuts along the path as lightweight
  `CutRef` objects. Gives O(depth) candidate sets per target with zero additional
  max-flow cost. Used by RCI for multi-candidate extraction.

**References**:
- Gusfield, D. (1990). Very simple methods for all pairs network flow analysis. *SIAM J. Comput.*, 19(1), 143-155.
- Dinitz, Y. (1970). Algorithm for solution of a problem of maximum flow in a network with power estimation. *Soviet Math. Doklady*, 11, 1277-1280.

## Subtour Elimination (SEC)

File: `src/sep/sec_separator.cpp`

For each customer t with y_t > tol, query the Gomory-Hu tree min-cut between
depot and t. The SEC requires:

```
x(delta(S)) >= 2 * y_t
```

where S is the target-side partition (not containing the depot).

For **s-t path** instances (source != target), sets containing the path target
need only 1 cut crossing (the path enters and terminates), so the RHS becomes
`y_t` instead of `2 * y_t`.

**Sparse form selection**: The separator counts nonzeros in two equivalent forms
and emits whichever is sparser:

| Form | Constraint | Nonzeros |
|---|---|---|
| Cut form | 2*y_t - x(delta(S)) <= 0 | \|delta(S)\| + 1 |
| Inside form | x(E(S)) - sum_{j in S\{t}} y_j <= 0 | \|E(S)\| + \|S\| - 1 |

The inside form is derived via degree substitution: x(delta({i})) = 2*y_i implies
x(delta(S)) = 2*sum(y_i, i in S) - 2*x(E(S)).

SEC is also used as a **feasibility check** on incumbent solutions found by HiGHS
heuristics (feasibility pump, rounding, sub-MIPs). Any integer solution with
a violated SEC is rejected. This uses tight tolerance (1e-6) vs the looser
fractional tolerance (default 0.1).

### Sub-MIP feasibility and separation

HiGHS primal heuristics (RENS, RINS) solve sub-MIPs that inherit our callbacks
but run their own presolve, which remaps columns. Two mechanisms handle this:

1. **Feasibility check with column mapping**: In `addIncumbent()`, sub-MIP
   solutions are expanded from reduced to original column space via
   `postSolveStack.undoPrimal()` before the SEC feasibility check. Without this,
   ~20-42% of sub-MIP incumbents violate SEC and corrupt the upper bound.

2. **Root-node SEC separation**: At the sub-MIP root node (`num_nodes == 0`),
   the LP solution is expanded to original space, SEC separation runs, and
   resulting cuts are translated back to reduced space. An inverse column map
   (`orig_to_reduced`) maps cut indices; eliminated columns have their fixed
   values subtracted from the cut RHS. Only SEC runs in sub-MIPs -- other
   separators are skipped to minimize overhead.

Sub-MIP separation is controlled by `--submip_separation true/false` (default: true).

## Rounded Capacity Inequalities (RCI)

File: `src/sep/rci_separator.cpp`

For a set S of customers with total demand d(S), define:

- Q_r = d(S) mod Q  (skip if Q_r = 0)
- k = ceil(d(S) / Q)  (skip if k <= 1)

The RCI in >= form:

```
x(delta(S)) - (2/Q_r) * sum(q_i * y_i, i in S) >= 2 * (k - d(S)/Q_r)
```

**Inside form** (via degree substitution, chosen when sparser):

```
x(E(S)) - sum_{i in S} (1 - q_i/Q_r) * y_i <= d(S)/Q_r - k
```

### Multi-candidate extraction

Unlike SEC which uses one min-cut per target, RCI uses `all_cuts_on_path(target)`
to extract **all** Gomory-Hu tree cuts along each target's path to root. This
yields O(depth) candidate sets per target with no extra max-flow calls. Duplicate
candidates (same GH tree node) are skipped via an `evaluated` bitmap.

### Add/drop local search

Each candidate set is refined by alternating add and drop phases:

1. **Drop phase**: Try removing each node from S; accept the best-improving removal.
   Repeat until no improvement.
2. **Add phase**: Try adding each node adjacent (in the support graph) to S;
   accept the best-improving addition. Repeat until no improvement.
3. Alternate drop/add until no improvement in either phase.

Candidates are evaluated and refined in parallel via TBB.

### Output cap

Results are sorted by violation (descending) and capped at max(5, n/5) cuts
per separator round (internal cap, before the per-separator top-k limit in
HiGHSBridge).

## Multistar / GLM Inequalities

File: `src/sep/multistar_separator.cpp`

The GLM inequality (Jepsen et al. eq 17) for a set S:

```
sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q) * x_e - (2/Q) * sum_{i in S} q_i * y_i >= 0
```

where t(e) is the endpoint of crossing edge e that lies **outside** S.

In <= form (for the HiGHS cut pool):

```
(2/Q) * sum_{i in S} q_i * y_i - sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q) * x_e <= 0
```

Uses one Gomory-Hu tree min-cut per target node, same as SEC.

## Rounded GLM (RGLM)

File: `src/sep/rglm_separator.h`, `src/sep/rglm_separator.cpp`

Strengthens GLM (multistar) inequalities using ceiling-based rounding, analogous
to how RCI strengthens the basic capacity inequality. Based on eq (24) of
Jepsen et al. 2014 with a_i = d_i, b = Q. Enabled via `--enable_rglm true`
(disabled by default).

### Background

The GLM inequality (eq 17) for a set S of customers uses the vehicle capacity Q
as the denominator in its coefficients. When the total demand is not a multiple
of Q, the required number of vehicles k = ceil(alpha/Q) leaves a fractional
remainder r = alpha mod Q. RGLM exploits this by replacing Q with the smaller
value r in the denominators, producing tighter coefficients.

The paper (p. 85) notes that RGLM is not separated directly -- instead, sets S
found by other separators (multistar, capacity) are tested for RGLM violations.
We follow this approach by reusing the Gomory-Hu tree min-cut sets.

### Inequality

For S subset of N (customer nodes, no depot), define:

- alpha = sum_{i in S} d_i + 2 * sum_{j not in S} d_j
- r = alpha mod Q (skip if r = 0, no rounding benefit)
- k = ceil(alpha / Q) (skip if k <= 1)

The RGLM constraint in >= form:

```
sum_{e in delta(S)} x_e  >=  2k - 2*beta/r
```

where beta = sum_{i in S} d_i(1-y_i) + sum_{j not in S} d_j(2 - sum_{e in E(j,S)} x_e).

Rearranging to <= form for the cut pool:

```
(2/r) * sum_{i in S} d_i * y_i
  + sum_{e in delta(S)} (2*d_out(e)/r - 1) * x_e
  <= 2 * (alpha/r - k)
```

where d_out(e) is the demand of the endpoint of e outside S.

### Key difference from GLM

| | GLM | RGLM |
|---|---|---|
| x-coefficient | -(1 - 2*d_out/Q) | 2*d_out/r - 1 |
| y-coefficient | 2*d_i/Q | 2*d_i/r |
| RHS | 0 | 2*(alpha/r - k) |
| Denominator | Q (capacity) | r = alpha mod Q |

Since r <= Q, the RGLM coefficients are at least as strong as GLM whenever
r > 0 and k > 1.

### Algorithm

Follows the same structure as `MultistarSeparator::separate()`:

1. Precompute d_total = sum of all customer demands (constant across targets).
2. For each target node (non-depot, y > tol):
   - Get min-cut from Gomory-Hu tree: S = {u : !reachable[u]}
   - Compute alpha = 2*d_total - d_S, r = alpha mod Q, k = ceil(alpha/Q)
   - Skip if r <= tol (no rounding benefit) or k <= 1
   - Compute cut_flow and weighted star_sum over crossing edges
   - Check violation: (2k - 2*beta/r) - cut_flow > tol
   - Emit cut in <= form with appropriate coefficients

### Usage

CLI:
```bash
./build/rcspp-solve instance.sppcc --enable_rglm true
```

C++ API:
```cpp
rcspp::SolverOptions opts = {{"enable_rglm", "true"}};
auto result = model.solve(opts);
```

Python:
```python
model.solve({"enable_rglm": "true"})
```

### Benchmark

Tested on 24 A/B/E SPPRCLIB instances (60s time limit):

| Metric | Value |
|---|---|
| Geo-mean speedup (ON vs OFF) | 1.03x |
| Wins ON / OFF / Tie | 10 / 8 / 6 |
| Best speedup | B-n57-k7-20 1.86x, B-n45-k6-54 1.46x |
| Worst slowdown | B-n52-k7-15 0.57x, A-n54-k7-149 0.71x |

Performance is neutral on aggregate. Helps on some B-class instances where the
rounding remainder is significant relative to capacity, but adds overhead on
instances where GLM cuts are already sufficient.

Disabled by default pending further tuning (e.g., filtering by violation
magnitude, limiting to sets with large r/Q ratio).

## Comb Inequalities

File: `src/sep/comb_separator.cpp`

Implements a BFS-based heuristic for comb separation (Jepsen et al. eq 13):

1. **Handle construction**: BFS from depot, adding nodes with x_e > threshold
   (tries thresholds 0.5 and 0.3).
2. **Tooth identification**: Crossing edges (one endpoint in handle, one outside)
   become candidate teeth. Each tooth contributes x_e - y_{outside}.
3. **Greedy selection**: Sort teeth by contribution, greedily select with no
   shared inside/outside nodes. Require >= 3 teeth, odd count.
4. **Violation check**: x(E(H)) + sum(x_{tooth}) - sum_{j in H}(y_j) - sum(y_{outside}) > (t-1)/2.

**Status**: Currently active but produces very few cuts in practice. The
formulation requires y-variable corrections not present in the standard CVRP
comb form.

## SeparationOracle (solver-independent)

File: `src/sep/separation_oracle.h`, `src/sep/separation_oracle.cpp`

The `SeparationOracle` class provides a solver-independent interface to all
separators. It handles support graph construction, Gomory-Hu tree computation,
and parallel separator execution — everything needed to generate cuts from an
LP relaxation solution.

### Usage

```cpp
#include "sep/separation_oracle.h"

rcspp::sep::SeparationOracle oracle(prob);
oracle.add_default_separators();  // SEC + RCI + Multistar + Comb

// In your solver's cut callback:
auto cuts = oracle.separate(x_values, y_values, x_offset, y_offset);
```

### Pipeline (per `separate()` call)

1. Build fractional support graph from edge LP values
2. Build Gomory-Hu tree (n-1 max-flow calls via Dinitz)
3. Build `SeparationContext` (LP values, problem, GH tree, offsets)
4. Run all registered separators **in parallel** (TBB task group)
5. For each separator: sort cuts by violation (descending), keep top-k
6. Merge and sort all cuts by violation

### Feasibility check

`oracle.is_feasible(x, y, x_off, y_off)` creates a temporary SEC separator
with tight tolerance (1e-6) and checks for violated subtour elimination
constraints. Used as a lazy-constraint callback on integer-feasible solutions.

### Parameters

| Method | Default | Description |
|---|---|---|
| `set_max_cuts_per_separator(k)` | 3 | Top-k most-violated cuts per separator per round (0 = unlimited) |
| `separate(..., tol)` | 0.1 | Violation tolerance for cut separation |

### Adding separators

Default set via `add_default_separators()`: SEC, RCI, Multistar, Comb.

Individual separators:
```cpp
oracle.add_separator(std::make_unique<rcspp::sep::SECSeparator>());
oracle.add_separator(std::make_unique<rcspp::sep::RCISeparator>());
oracle.add_separator(std::make_unique<rcspp::sep::MultistarSeparator>());
oracle.add_separator(std::make_unique<rcspp::sep::CombSeparator>());
oracle.add_separator(std::make_unique<rcspp::sep::RGLMSeparator>());
```

### Cut format

Each `Cut` returned by `separate()`:

| Field | Type | Description |
|---|---|---|
| `indices` | `vector<int32_t>` | Column indices in the LP |
| `values` | `vector<double>` | Coefficients |
| `rhs` | `double` | Right-hand side (<= form) |
| `violation` | `double` | Violation amount |

## Cut Management (HiGHS)

Managed by `HiGHSBridge::install_separators()` (`src/model/highs_bridge.cpp`).
The HiGHS integration uses `SeparationOracle` internally for feasibility checking
and wraps the same separation pipeline for HiGHS's user separator callback.

### Per-round pipeline

1. Build fractional support graph from LP relaxation
2. Build Gomory-Hu tree (shared, n-1 max-flows)
3. Run all separators **in parallel** (TBB task group)
4. For each separator: sort cuts by violation (descending), keep top-k

### Parameters

| Parameter | CLI flag | Default | Description |
|---|---|---|---|
| max_cuts_per_sep | `--max_cuts_per_separator N` | 3 | Top-k most-violated cuts per separator per round |
| separation_tol | `--separation_tol X` | 0.1 | Fractional violation threshold |
| separation_interval | `--separation_interval N` | 1 | Run separation every N-th callback (amortization) |
| enable_rglm | `--enable_rglm true` | false | Enable RGLM separator |
| submip_separation | `--submip_separation true` | true | SEC separation at sub-MIP root node |

### Tolerances

Two separate tolerances:

- **Fractional** (frac_tol_, default 0.1): Used during LP separation. Higher
  threshold avoids weak cuts that slow the LP. Jepsen et al. (2008) found 0.4
  optimal with CPLEX; HiGHS needs lower since it has fewer built-in cuts.
- **Integral** (int_tol_, 1e-6): Used for feasibility checking of incumbent
  solutions. Tight to avoid rejecting valid solutions due to float noise.

## References

- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). A branch-and-cut algorithm for the capacitated profitable tour problem. *Discrete Optimization*, 14, 78-96.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2008). Subset-row inequalities applied to the vehicle-routing problem with time windows. *Operations Research*, 56(2), 497-511.
