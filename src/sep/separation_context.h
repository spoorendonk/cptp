#pragma once

#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace cptp {
class Problem;
class digraph;
class gomory_hu_tree;
}  // namespace cptp

namespace cptp::sep {

class LargestSideCache;

/// Default fractional separation tolerance.
/// Jepsen et al. (2014) found 0.4 optimal with CPLEX's built-in cuts.
/// HiGHS needs more custom cuts, so we use a lower threshold.
inline constexpr double kDefaultFracTol = 0.1;

/// Solver-independent context passed to each separator.
struct SeparationContext {
    const Problem& problem;
    std::span<const double> x_values;  // edge variables
    std::span<const double> y_values;  // node variables
    int32_t x_offset;                  // LP column offset for edge vars
    int32_t y_offset;                  // LP column offset for node vars
    double tol = 1e-6;                 // cut violation tolerance
    const gomory_hu_tree* flow_tree = nullptr;

    /// The support graph `flow_tree` was built on (arcs both ways, x_e as
    /// capacity), for separators that run their own max-flows on it, such as
    /// the multistar separator's largest-side option (#168). Null when the
    /// caller does not provide it.
    const digraph* support_graph = nullptr;
    const std::vector<double>* support_capacity = nullptr;

    /// The largest sides on `support_graph`, shared by the separators that
    /// use them (#175). Null when the caller does not provide it: each such
    /// separator then runs its own max-flows on the support graph.
    LargestSideCache* largest_sides = nullptr;

    /// Current incumbent upper bound (for cost-based separators).
    double upper_bound = std::numeric_limits<double>::infinity();

    /// All-pairs shortest path bounds: flat n×n matrix, d(s,v) =
    /// all_pairs[s*n+v]. Empty when not available. Used by SPISeparator.
    /// `= {}` keeps gcc's -Wmissing-field-initializers quiet for the
    /// designated initializers that leave it out.
    std::span<const double> all_pairs = {};  // NOLINT(readability-redundant-member-init)

    /// Threads a separator may use for its own loop, the calling thread
    /// included; 1 runs on the calling thread. The bridge sets it to what the
    /// solve's `threads` budget leaves after the threads running the
    /// separators, and 1 under HiGHS's parallel tree search, whose workers
    /// already use the budget. Every separator sees the same value, so the
    /// budget holds only while a single separator (RCI) reads it (#159).
    int32_t threads = 1;
};

}  // namespace cptp::sep
