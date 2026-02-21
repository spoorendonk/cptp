#pragma once

#include <cstdint>
#include <span>

namespace cptp {
class Problem;
class gomory_hu_tree;
}

namespace cptp::sep {

/// Default fractional separation tolerance.
/// Jepsen et al. (2008) found 0.4 optimal with CPLEX's built-in cuts.
/// HiGHS needs more custom cuts, so we use a lower threshold.
inline constexpr double kDefaultFracTol = 0.1;

/// Solver-independent context passed to each separator.
struct SeparationContext {
    const Problem& problem;
    std::span<const double> x_values;  // edge variables
    std::span<const double> y_values;  // node variables
    int32_t x_offset;  // LP column offset for edge vars
    int32_t y_offset;  // LP column offset for node vars
    double tol = 1e-6;                 // cut violation tolerance
    const gomory_hu_tree* flow_tree = nullptr;
};

}  // namespace cptp::sep
