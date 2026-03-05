#pragma once

#include <limits>
#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "Highs.h"
#include "core/problem.h"
#include "core/solution.h"
#include "sep/separation_context.h"
#include "sep/separation_oracle.h"
#include "sep/separator.h"
#include "util/logger.h"
namespace cptp {

/// HiGHS options as string key-value pairs, forwarded directly.
using SolverOptions = std::vector<std::pair<std::string, std::string>>;

/// Strategy for Lagrangian reduced-cost fixing in the propagator.
enum class RCFixingStrategy { off, root_only, on_ub_improvement, periodic, adaptive };

struct RCFixingSettings {
    RCFixingStrategy strategy = RCFixingStrategy::off;
    int32_t periodic_interval = 100;
    bool fix_to_one = false;
};

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
    int32_t num_nodes() const { return num_nodes_; }

    /// Set separation skip interval: separate every N-th callback invocation.
    /// 1 = every round (default), 2 = every other round, etc.
    void set_separation_interval(int32_t interval) { separation_interval_ = interval; }

    /// Max cuts to add per separator per round (0 = unlimited).
    /// Jepsen et al. recommend 1: add only the most-violated cut.
    void set_max_cuts_per_separator(int32_t max_cuts) { max_cuts_per_sep_ = max_cuts; }
    /// Global cut cap: pool all cuts, sort by violation, add top N (0=unlimited).
    void set_max_cuts_per_round(int32_t max_cuts) { max_cuts_per_round_ = max_cuts; }
    void set_separator_max_cuts(const std::string& name, int32_t max_cuts) {
        per_separator_max_cuts_[name] = max_cuts;
    }
    void set_separator_min_violation(const std::string& name, double min_violation) {
        per_separator_min_violation_[name] = min_violation;
    }

    /// Enable SEC separation inside HiGHS sub-MIPs (RENS/RINS).
    /// Requires column mapping via undoPrimal to translate between reduced/original space.
    void set_submip_separation(bool enable) { submip_separation_ = enable; }

    /// Set upper bound for edge elimination preprocessing.
    void set_upper_bound(double ub) { upper_bound_ = ub; }

    /// Configure Lagrangian reduced-cost fixing in the propagator.
    void set_rc_fixing(RCFixingSettings settings) { rc_settings_ = settings; }
    void set_edge_elimination(bool enable) { edge_elimination_enabled_ = enable; }
    void set_edge_elimination_nodes(bool enable) { edge_elimination_nodes_ = enable; }

    /// Set pre-computed labeling bounds for edge elimination and propagation.
    /// f = forward bounds (source → v), b = backward bounds (v → target).
    /// correction = profit(source) if source==target (tour), 0 if source!=target (path).
    void set_labeling_bounds(std::vector<double> f, std::vector<double> b, double correction);

    /// Set all-pairs labeling bounds for stronger domain propagation (Trigger B).
    /// dist is a flat n×n matrix: d(s,v) = dist[s*n + v].
    void set_all_pairs_bounds(std::vector<double> dist);

    void set_interrupt_flag(std::shared_ptr<std::atomic<bool>> flag) {
        interrupt_flag_ = std::move(flag);
    }

    /// Install domain propagator that fixes edges during B&C based on labeling bounds.
    void install_propagator();

    /// Install LP-guided primal heuristic callback (kCallbackMipUserSolution).
    void install_heuristic_callback();

    /// Enable/disable heuristic callback.
    void set_heuristic_callback(bool enable) { heuristic_callback_ = enable; }

    /// Set time budget for each heuristic callback invocation (ms).
    void set_heuristic_budget_ms(double ms) { heuristic_budget_ms_ = ms; }

    /// Set heuristic strategy: 0=all, 1=LP-threshold, 2=RINS, 3=neighborhood.
    void set_heuristic_strategy(int s) { heuristic_strategy_ = s; }
    void set_heuristic_node_interval(int64_t interval) { heuristic_node_interval_ = interval; }
    void set_heuristic_deterministic_restarts(int32_t restarts) {
        heuristic_deterministic_restarts_ = restarts;
    }
    /// Control bounds-based propagation (Trigger A/B) independently of RC fixing.
    void set_bounds_propagation(bool enable) { bounds_propagation_ = enable; }

    /// LP-guided heuristic threshold parameters.
    void set_heu_lpg_edge_threshold(double v) { heu_lpg_edge_threshold_ = v; }
    void set_heu_lpg_node_threshold(double v) { heu_lpg_node_threshold_ = v; }
    void set_heu_lpg_lp_threshold(double v) { heu_lpg_lp_threshold_ = v; }
    void set_heu_lpg_seed_threshold(double v) { heu_lpg_seed_threshold_ = v; }

    void set_fixed_y(std::vector<int8_t> fixed_y) { fixed_y_ = std::move(fixed_y); }

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
    int32_t max_cuts_per_round_ = 0;   // global cut cap per round (0 = unlimited)
    std::unordered_map<std::string, int32_t> per_separator_max_cuts_;
    std::unordered_map<std::string, double> per_separator_min_violation_;
    bool submip_separation_ = true;    // SEC separation at sub-MIP root node
    double upper_bound_ = std::numeric_limits<double>::infinity();
    bool edge_elimination_enabled_ = true;
    bool edge_elimination_nodes_ = true;
    bool bounds_propagation_ = true;

    // LP-guided heuristic thresholds
    double heu_lpg_edge_threshold_ = 0.1;
    double heu_lpg_node_threshold_ = 0.5;
    double heu_lpg_lp_threshold_ = 0.1;
    double heu_lpg_seed_threshold_ = 0.3;

    // Labeling bounds for edge elimination and propagation
    std::vector<double> fwd_bounds_;   // f[v]: forward bounds (source → v)
    std::vector<double> bwd_bounds_;   // b[v]: backward bounds (v → target)
    double correction_ = 0.0;         // profit(source) if tour, 0 if path
    std::vector<int8_t> fixed_y_;      // optional fixed y-vars: -1 unset, 0/1 fixed

    // All-pairs labeling bounds for stronger Trigger B propagation
    std::vector<double> all_pairs_;    // flat n×n: d(s,v) = all_pairs_[s*n+v]
    std::shared_ptr<std::atomic<bool>> interrupt_flag_;

    // RC fixing settings and statistics
    RCFixingSettings rc_settings_;
    std::shared_ptr<int64_t> rc_fix0_count_;
    std::shared_ptr<int64_t> rc_fix1_count_;
    std::shared_ptr<int64_t> rc_label_runs_;
    std::shared_ptr<int64_t> rc_callback_runs_;
    std::shared_ptr<double> rc_time_seconds_;

    // Propagator statistics (shared with callback lambda)
    std::shared_ptr<int64_t> propagator_fixings_;
    std::shared_ptr<int64_t> sweep_fixings_;
    std::shared_ptr<int64_t> chain_fixings_;
    std::shared_ptr<int64_t> ub_improvements_;
    std::shared_ptr<int64_t> propagator_calls_;
    std::shared_ptr<double> propagator_time_seconds_;

    // Cut statistics (updated from callback)
    mutable std::map<std::string, SeparatorStats> separator_stats_;
    mutable int64_t total_cuts_ = 0;
    mutable int64_t separation_rounds_ = 0;
    mutable int64_t separation_calls_ = 0;  // total callback invocations

    // Heuristic callback: LP-guided primal heuristic
    bool heuristic_callback_ = true;
    double heuristic_budget_ms_ = 20.0;
    int heuristic_strategy_ = 0;  // 0=all, 1=threshold, 2=RINS, 3=neighborhood
    int64_t heuristic_node_interval_ = 200;
    int32_t heuristic_deterministic_restarts_ = 32;

    // Cached LP relaxation from separator callback (shared with heuristic)
    std::shared_ptr<std::vector<double>> cached_x_lp_;
    std::shared_ptr<std::vector<double>> cached_y_lp_;
    std::shared_ptr<std::mutex> lp_cache_mutex_;

    // Heuristic callback statistics
    std::shared_ptr<int64_t> heuristic_calls_;
    std::shared_ptr<int64_t> heuristic_solutions_;
    std::shared_ptr<int64_t> heuristic_work_units_;
    std::shared_ptr<double> heuristic_time_seconds_;
};

}  // namespace cptp
