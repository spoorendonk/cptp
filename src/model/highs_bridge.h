#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Highs.h"
#include "core/problem.h"
#include "core/solution.h"
#include "sep/separation_context.h"
#include "sep/separator.h"

namespace cptp {

/// HiGHS options as string key-value pairs, forwarded directly.
using SolverOptions = std::vector<std::pair<std::string, std::string>>;

/// Wires the CPTP formulation and custom separators into HiGHS.
class HiGHSBridge {
 public:
    HiGHSBridge(const Problem& prob, Highs& highs,
                 double frac_tol = sep::kDefaultFracTol);
    ~HiGHSBridge();

    void add_separator(std::unique_ptr<sep::Separator> sep);

    /// Build the IP formulation: variables + core constraints.
    void build_formulation();

    /// Register separators into HiGHS's cutpool via HighsUserSeparator callback.
    void install_separators();

    /// Extract solution after solve.
    SolveResult extract_result() const;

    int32_t x_offset() const { return 0; }
    int32_t y_offset() const { return num_edges_; }

    /// Set separation skip interval: separate every N-th callback invocation.
    /// 1 = every round (default), 2 = every other round, etc.
    void set_separation_interval(int32_t interval) { separation_interval_ = interval; }

    /// Max cuts to add per separator per round (0 = unlimited).
    /// Jepsen et al. recommend 1: add only the most-violated cut.
    void set_max_cuts_per_separator(int32_t max_cuts) { max_cuts_per_sep_ = max_cuts; }

 private:
    /// Order the tour nodes by following edges from depot.
    std::vector<int32_t> order_tour(const std::vector<int32_t>& visited_nodes,
                                     const std::vector<int32_t>& active_edges) const;

    const Problem& prob_;
    Highs& highs_;
    int32_t num_edges_;
    int32_t num_nodes_;
    double frac_tol_;   // violation tolerance for fractional separation
    double int_tol_;    // violation tolerance for integral feasibility check
    std::vector<std::unique_ptr<sep::Separator>> separators_;

    // Amortized separation: skip rounds to reduce overhead
    int32_t separation_interval_ = 1;  // 1 = every round (no skipping)
    int32_t max_cuts_per_sep_ = 3;     // max cuts per separator per round (0 = unlimited)

    // Cut statistics (updated from callback)
    mutable std::map<std::string, SeparatorStats> separator_stats_;
    mutable int64_t total_cuts_ = 0;
    mutable int64_t separation_rounds_ = 0;
    mutable int64_t separation_calls_ = 0;  // total callback invocations
};

}  // namespace cptp
