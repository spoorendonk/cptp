#pragma once

#include "sep/separator.h"

#include <cstdint>
#include <span>
#include <vector>

namespace cptp {
class Problem;
}  // namespace cptp

namespace cptp::sep {

/// Lift the cover inequality sum_{i in C} y_i <= |C| - 1 against every
/// bounded capacity row (#166, #146 Q2): the nodes of `order` are lifted one
/// after another, alpha_j = (|C| - 1) - U_j, where U_j is the least over the
/// rows of the exact single-row lifting value given the coefficients lifted so
/// far (each an upper bound on the all-row lifting problem, so the result stays
/// valid). A node that cannot be served alone gets |C| - 1. Returns alpha per
/// node: 1 on the cover, the lifted value on `order`, 0 elsewhere. The cover
/// needs |C| >= 2 and must be a cover of some row; the depot must not be in
/// either span.
std::vector<int32_t> lift_cover_all_rows(const Problem& prob, std::span<const int32_t> cover,
                                         std::span<const int32_t> order);

/// Cross-row lifted cover w-GLM separator (#166, theory on #146 Q2-Q4 §1-§2).
/// For each candidate set S (each target's Gomory-Hu min-cut side, its
/// largest min-cut side when the context has the support graph, and V minus
/// the depot), with z_i = y_i in S and z_j = x(j : S) outside, it finds per
/// bounded row a cover over z, lifts it against all rows, and emits the
/// w-GLM
///   x(delta(S)) >= 2 sum_{i in S} w_i y_i + 2 sum_{e in delta(S)} w_{t(e)} x_e
/// with w = alpha / (|C| - 1), w_0 = 0, when the LP point violates it.
class XCoverSeparator : public Separator {
public:
    [[nodiscard]] std::string name() const override { return "XCover"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace cptp::sep
