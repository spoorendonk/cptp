#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace cptp::sep {

/// A coefficient is near zero when |a_j| <= kNearZeroCoefficient * max_k |a_k|.
/// Our separators' coefficients are 1, 2, or ratios of demands to
/// capacities (RCI's 2*q_i/Q, a GLM's 1 - 2*q_i/Q), so the smallest one a
/// separator means is about min q_i / Q, far above 1e-9 on any instance we
/// solve; floating-point noise in such a ratio (a GLM's 1 - 2*q/Q for q = Q/2
/// can come out as ~1e-17 instead of 0) is about 1e-16 relative, far below
/// it. HiGHS's own cut cleanup (postprocessCut) drops below 1e-4 relative,
/// which would also drop small RCI coefficients on instances with a large Q.
inline constexpr double kNearZeroCoefficient = 1e-9;

/// Drops the near-zero coefficients of the cut sum_k val[k] * x[idx[k]] <= rhs
/// (#191) before it goes to HiGHS, whose clique extraction divides by every
/// coefficient: a zero gives a NaN bound (#188), a tiny one a huge implied
/// bound. Each dropped term moves to the rhs at the bound that relaxes it,
/// so the cut stays valid: rhs -= a_j * lower[j] for a_j > 0 and
/// rhs -= a_j * upper[j] for a_j < 0. The bounds must be valid for the whole
/// tree. Exact zeros drop without touching the rhs.
///
/// Returns false when the cut should not be added: a term to drop whose
/// relaxing bound is infinite, nothing left, or, when a nonzero term was
/// dropped, a violation at `lp_point` no longer above `min_violation`.
/// `idx`, `val` and `rhs` are then unspecified.
[[nodiscard]] bool drop_near_zero_coefficients(std::vector<int32_t>& idx, std::vector<double>& val,
                                               double& rhs, std::span<const double> lower,
                                               std::span<const double> upper,
                                               std::span<const double> lp_point,
                                               double min_violation);

}  // namespace cptp::sep
