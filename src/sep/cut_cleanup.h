#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace cptp::sep {

/// A coefficient is near zero when |a_j| <= kNearZeroCoefficient, absolute.
/// Our separators' coefficients are 1, 2, or ratios of demands to
/// capacities (RCI's 2*q_i/Q_r and 1 - q_i/Q_r with Q_r = d(S) mod Q, a
/// GLM's 1 - 2*q_i/Q), so the smallest one a separator means is about
/// min q_i / Q or 1 / Q_r, far above 1e-9 on any instance we solve, while
/// floating-point noise in such an O(1) expression (a GLM's 1 - 2*q/Q for
/// q = Q/2 can come out as ~1e-17 instead of 0) is about 1e-16, far below
/// it. Absolute rather than relative to max|a|: RCI's 2*q_i/Q_r grows
/// without bound as Q_r gets small, and a relative threshold would then drop
/// the cut's +-1 edge terms, losing violation of order 1. HiGHS's own cut
/// cleanup (postprocessCut), at 1e-4 relative, would also drop small RCI
/// coefficients on instances with a large Q.
inline constexpr double kNearZeroCoefficient = 1e-9;

/// Drops the near-zero coefficients of the cut sum_k val[k] * x[idx[k]] <= rhs
/// (#191) before it goes to HiGHS, whose clique extraction divides by every
/// coefficient: a zero gives a NaN bound (#188), a tiny one a huge implied
/// bound. Each dropped term moves to the rhs at the bound that relaxes it,
/// so the cut stays valid: rhs -= a_j * lower[j] for a_j > 0 and
/// rhs -= a_j * upper[j] for a_j < 0. The bounds must be valid for the whole
/// tree. Exact zeros drop without touching the rhs.
///
/// Returns how much less the cut is violated at `lp_point` (activity minus
/// rhs) after the drop, sum_j a_j * (lp_point[j] - bound_j) over the dropped
/// terms, which is >= 0 for a point within the bounds; the caller decides
/// whether the cut still cuts off the point by enough. Returns nullopt when
/// the cut should not be added: a coefficient that is not finite, a term to
/// drop whose relaxing bound is infinite, or nothing left. `idx`, `val` and
/// `rhs` are then unspecified.
[[nodiscard]] std::optional<double> drop_near_zero_coefficients(
    std::vector<int32_t>& idx, std::vector<double>& val, double& rhs, std::span<const double> lower,
    std::span<const double> upper, std::span<const double> lp_point);

}  // namespace cptp::sep
