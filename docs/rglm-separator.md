# RGLM (Rounded GLM) Separator

File: `src/sep/rglm_separator.h`, `src/sep/rglm_separator.cpp`

## Purpose

Strengthens GLM (multistar) inequalities using ceiling-based rounding, analogous
to how RCI (eq 21) strengthens the basic capacity inequality. Based on eq (24)
of Jepsen et al. 2014 with a_i = d_i, b = Q.

Opt-in via `--enable_rglm true` (disabled by default).

## Background

The GLM inequality (eq 17) for a set S of customers uses the vehicle capacity Q
as the denominator in its coefficients. When the total demand is not a multiple
of Q, the required number of vehicles k = ceil(alpha/Q) leaves a fractional
remainder r = alpha mod Q. RGLM exploits this by replacing Q with the smaller
value r in the denominators, producing tighter coefficients.

The paper (p. 85) notes that RGLM is not separated directly -- instead, sets S
found by other separators (multistar, capacity) are tested for RGLM violations.
We follow this approach by reusing the Gomory-Hu tree min-cut sets.

## Inequality

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

## Key difference from GLM

| | GLM | RGLM |
|---|---|---|
| x-coefficient | -(1 - 2*d_out/Q) | 2*d_out/r - 1 |
| y-coefficient | 2*d_i/Q | 2*d_i/r |
| RHS | 0 | 2*(alpha/r - k) |
| Denominator | Q (capacity) | r = alpha mod Q |

Since r <= Q, the RGLM coefficients are at least as strong as GLM whenever
r > 0 and k > 1.

## Algorithm

Follows the same structure as `MultistarSeparator::separate()`:

1. Precompute d_total = sum of all customer demands (constant across targets).
2. For each target node (non-depot, y > tol):
   - Get min-cut from Gomory-Hu tree: S = {u : !reachable[u]}
   - Compute alpha = 2*d_total - d_S, r = alpha mod Q, k = ceil(alpha/Q)
   - Skip if r <= tol (no rounding benefit) or k <= 1
   - Compute cut_flow and weighted star_sum over crossing edges
   - Check violation: (2k - 2*beta/r) - cut_flow > tol
   - Emit cut in <= form with appropriate coefficients

## Usage

CLI:
```bash
./build/cptp-solve instance.sppcc --enable_rglm true
```

C++ API:
```cpp
cptp::SolverOptions opts = {{"enable_rglm", "true"}};
auto result = model.solve(opts);
```

Python:
```python
model.solve({"enable_rglm": "true"})
```

## Benchmark

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

## References

- Jepsen et al. 2014, "A branch-and-cut algorithm for the capacitated profitable
  tour problem", Discrete Optimization 14, pp. 78-96. Eq (24), Lemma 1.
