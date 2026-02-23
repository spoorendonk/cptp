#pragma once

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Highs.h"
#include "core/problem.h"
#include "core/solution.h"
#include "sep/separation_context.h"
#include "sep/separation_oracle.h"
#include "sep/separator.h"
#include "util/logger.h"

namespace rcspp {

/// HiGHS options as string key-value pairs, forwarded directly.
using SolverOptions = std::vector<std::pair<std::string, std::string>>;

/// Wires the CPTP formulation and custom separators into HiGHS.
class HiGHSBridge {
 public:
    HiGHSBridge(const Problem& prob, Highs& highs, Logger& logger,
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

    /// Enable SEC separation inside HiGHS sub-MIPs (RENS/RINS).
    /// Requires column mapping via undoPrimal to translate between reduced/original space.
    void set_submip_separation(bool enable) { submip_separation_ = enable; }

    /// Set upper bound for edge elimination preprocessing.
    void set_upper_bound(double ub) { upper_bound_ = ub; }

    /// Set pre-computed labeling bounds for edge elimination and propagation.
    /// f = forward bounds (source → v), b = backward bounds (v → target).
    /// correction = profit(source) if source==target (tour), 0 if source!=target (path).
    void set_labeling_bounds(std::vector<double> f, std::vector<double> b, double correction);

    /// Set all-pairs labeling bounds for stronger domain propagation (Trigger B).
    /// dist is a flat n×n matrix: d(s,v) = dist[s*n + v].
    void set_all_pairs_bounds(std::vector<double> dist);

    /// Install domain propagator that fixes edges during B&C based on labeling bounds.
    void install_propagator();

 private:
    /// Order the visited nodes by following edges from source.
    std::vector<int32_t> order_path(const std::vector<int32_t>& visited_nodes,
                                     const std::vector<int32_t>& active_edges) const;

    const Problem& prob_;
    Highs& highs_;
    Logger& logger_;
    int32_t num_edges_;
    int32_t num_nodes_;
    double frac_tol_;   // violation tolerance for fractional separation
    double int_tol_;    // violation tolerance for integral feasibility check
    std::vector<std::unique_ptr<sep::Separator>> separators_;

    // Amortized separation: skip rounds to reduce overhead
    int32_t separation_interval_ = 1;  // 1 = every round (no skipping)
    int32_t max_cuts_per_sep_ = 3;     // max cuts per separator per round (0 = unlimited)
    bool submip_separation_ = true;    // SEC separation at sub-MIP root node
    double upper_bound_ = std::numeric_limits<double>::infinity();

    // Labeling bounds for edge elimination and propagation
    std::vector<double> fwd_bounds_;   // f[v]: forward bounds (source → v)
    std::vector<double> bwd_bounds_;   // b[v]: backward bounds (v → target)
    double correction_ = 0.0;         // profit(source) if tour, 0 if path

    // All-pairs labeling bounds for stronger Trigger B propagation
    std::vector<double> all_pairs_;    // flat n×n: d(s,v) = all_pairs_[s*n+v]

    // Propagator statistics (shared with callback lambda)
    std::shared_ptr<int64_t> propagator_fixings_;
    std::shared_ptr<int64_t> sweep_fixings_;
    std::shared_ptr<int64_t> chain_fixings_;
    std::shared_ptr<int64_t> ub_improvements_;
    std::shared_ptr<int64_t> propagator_calls_;

    // Cut statistics (updated from callback)
    mutable std::map<std::string, SeparatorStats> separator_stats_;
    mutable int64_t total_cuts_ = 0;
    mutable int64_t separation_rounds_ = 0;
    mutable int64_t separation_calls_ = 0;  // total callback invocations
};

}  // namespace rcspp
