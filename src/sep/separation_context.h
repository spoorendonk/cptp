#pragma once

#include <cstdint>
#include <limits>
#include <span>

namespace rcspp {
class Problem;
class gomory_hu_tree;
}

namespace rcspp::sep {

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

    /// Current incumbent upper bound (for cost-based separators).
    double upper_bound = std::numeric_limits<double>::infinity();

    /// All-pairs shortest path bounds: flat n×n matrix, d(s,v) = all_pairs[s*n+v].
    /// Empty when not available. Used by SPISeparator.
    std::span<const double> all_pairs = {};
};

}  // namespace rcspp::sep
