#pragma once

#include "core/problem.h"
#include "core/solution.h"
#include "Highs.h"
#include "model/hyperplane.h"
#include "preprocess/edge_elimination.h"
#include "sep/separation_context.h"
#include "sep/separation_oracle.h"
#include "sep/separator.h"
#include "util/logger.h"

#include <algorithm>
#include <atomic>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
namespace cptp {

/// HiGHS options as string key-value pairs, forwarded directly.
using SolverOptions = std::vector<std::pair<std::string, std::string>>;

/// Strategy for Lagrangian reduced-cost fixing in the propagator.
enum class RCFixingStrategy { off, root_only, on_ub_improvement, periodic, adaptive };

struct RCFixingSettings {
    RCFixingStrategy strategy = RCFixingStrategy::adaptive;
    int32_t periodic_interval = 100;
    bool fix_to_one = false;
};

/// Aggregated statistics from HiGHS callback lambdas.
/// Single allocation shared between the bridge and its callback closures.
struct CallbackStats {
    // Propagator
    int64_t propagator_fixings = 0;
    int64_t sweep_fixings = 0;
    int64_t chain_fixings = 0;
    int64_t sweep_node_fixings = 0;
    int64_t chain_node_fixings = 0;
    int64_t ub_improvements = 0;
    int64_t propagator_calls = 0;
    double propagator_time_seconds = 0.0;

    // RC fixing
    int64_t rc_fix0_count = 0;
    int64_t rc_fix1_count = 0;
    int64_t rc_label_runs = 0;
    int64_t rc_callback_runs = 0;
    double rc_time_seconds = 0.0;

    // Heuristic
    int64_t heuristic_calls = 0;
    int64_t heuristic_solutions = 0;
    double heuristic_time_seconds = 0.0;
};

/// The propagator's mutable state for one HiGHS tree worker, keyed by the
/// worker's LP relaxation and the solve's restart count: a worker runs its
/// nodes one after another, so its state needs no lock, while workers of a
/// parallel tree search run at once. A restart recreates the workers at the
/// same addresses, and the restart count keeps them from inheriting state.
/// This relies on HiGHS v1.16 recreating workers only at a restart; recheck it
/// on a HiGHS bump.
struct PropagatorWorkerState {
    double last_ub = std::numeric_limits<double>::infinity();
    std::vector<bool> processed_edge;
    CallbackStats stats;  // propagator and RC fields only, ub_improvements not
};

struct PropagatorStates {
    std::mutex mutex;  // guards by_lp, not the states
    // (&lp, mipdata_->numRestarts) -> that worker's state.
    std::map<std::pair<const void*, int64_t>, std::unique_ptr<PropagatorWorkerState>> by_lp;
    // What the RC fixing strategies decide on, solve-wide so that it neither
    // multiplies with the workers nor resets at a restart: the best UB seen by
    // any worker and how often it improved, the propagator calls and RC runs,
    // and the adaptive strategy's low-yield streak.
    std::atomic<double> best_ub = std::numeric_limits<double>::infinity();
    std::atomic<int64_t> ub_improvements = 0;
    std::atomic<int64_t> propagator_calls = 0;
    std::atomic<int64_t> rc_callback_runs = 0;
    std::atomic<int32_t> rc_low_yield_streak = 0;
    std::atomic<bool> rc_adaptive_disabled = false;
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

    /// The solve's thread budget (the resolved `threads` option, >= 1). Work
    /// in the callbacks (separation, reduced-cost fixing, the LP-guided
    /// heuristic) uses at most this many threads, the calling one included,
    /// and only the calling one under HiGHS's parallel tree search (#159).
    void set_threads(int32_t threads) { threads_ = std::max<int32_t>(1, threads); }

    /// Max cuts to add per separator per round (0 = unlimited).
    /// Jepsen et al. recommend 1: add only the most-violated cut.
    void set_max_cuts_per_separator(int32_t max_cuts) { max_cuts_per_sep_ = max_cuts; }
    /// Global cut cap: pool all cuts, sort by violation, add top N (0=unlimited).
    void set_max_cuts_per_round(int32_t max_cuts) { max_cuts_per_round_ = max_cuts; }
    void set_separator_max_cuts(const std::string& name, int32_t max_cuts) {
        per_separator_max_cuts_[name] = max_cuts;
    }

    /// Enable SEC separation inside HiGHS sub-MIPs (RENS/RINS).
    /// Requires column mapping via undoPrimal to translate between
    /// reduced/original space.
    void set_submip_separation(bool enable) { submip_separation_ = enable; }

    /// Set upper bound for edge elimination preprocessing.
    void set_upper_bound(double ub) { upper_bound_ = ub; }

    /// Warm-start solution injected through the user-solution callback (on the
    /// after-setup query) instead of HiGHS setSolution(). Routing it through the
    /// callback avoids a HiGHS 1.15 setSolution()+MipUserSolution interaction
    /// that corrupts the search and produces wrong proven optima.
    void set_warm_start_solution(std::vector<double> sol) { warm_start_solution_ = std::move(sol); }

    /// Configure Lagrangian reduced-cost fixing in the propagator.
    void set_rc_fixing(RCFixingSettings settings) { rc_settings_ = settings; }
    void set_labeling_max_queue_pops(int64_t limit) { labeling_max_queue_pops_ = limit; }
    void set_label_rows(preprocess::LabelRows rows) { label_rows_ = rows; }
    /// When the in-tree vector labeling stops: the solve's time limit (#115).
    void set_labeling_deadline(preprocess::LabelingDeadline deadline) {
        labeling_deadline_ = deadline;
    }
    void set_edge_elimination(bool enable) { edge_elimination_enabled_ = enable; }

    /// True if all objective coefficients are integer-valued (detected in
    /// build_formulation). Used to tighten LP bounds via ceiling in RC-fixing.
    bool obj_is_integer() const { return obj_is_integer_; }
    void set_edge_elimination_nodes(bool enable) { edge_elimination_nodes_ = enable; }

    /// Set pre-computed labeling bounds for edge elimination and propagation.
    /// f = forward bounds (source → v), b = backward bounds (v → target).
    /// correction = profit(source) if source==target (tour), 0 if source!=target
    /// (path).
    void set_labeling_bounds(std::vector<double> f, std::vector<double> b, double correction);

    /// Per-edge lower bounds (preprocess::edge_lower_bounds, #144) that edge
    /// elimination and the propagator's sweep use alongside f and b: an
    /// edge's bound is the larger of the two. Empty: f and b alone.
    void set_edge_lower_bounds(std::vector<double> edge_lb) { edge_lb_ = std::move(edge_lb); }

    /// Set all-pairs labeling bounds for stronger domain propagation (Trigger B).
    /// dist is a flat n×n matrix: d(s,v) = dist[s*n + v].
    void set_all_pairs_bounds(std::vector<double> dist);

    void set_interrupt_flag(std::shared_ptr<std::atomic<bool>> flag) {
        interrupt_flag_ = std::move(flag);
    }

    /// Install domain propagator that fixes edges during B&C based on labeling
    /// bounds.
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

    /// Hyperplanes to branch on, set before build_formulation: each becomes an
    /// integer column z_h after the x and y columns, linked by the row
    /// z_h - sum(a_j * y_j) = 0, and HiGHS branches on z_h as on any column
    /// (#137). The x and y offsets do not move.
    void set_hyperplanes(std::vector<detail::Hyperplane> hyperplanes) {
        hyperplanes_ = std::move(hyperplanes);
    }
    int32_t num_hyperplane_columns() const { return static_cast<int32_t>(hyperplanes_.size()); }

    /// Pads `sol` (x and y values, at least) to the x and y columns and appends
    /// each hyperplane column's z = a . y, so it is a full HiGHS solution.
    void complete_solution(std::vector<double>& sol) const;

private:
    /// Order the visited nodes by following edges from source.
    std::vector<int32_t> order_path(const std::vector<int32_t>& visited_nodes,
                                    const std::vector<int32_t>& active_edges) const;

    const Problem& prob_;
    Highs& highs_;
    Logger& logger_;
    int32_t num_edges_;
    int32_t num_nodes_;
    double frac_tol_;       // violation tolerance for fractional separation
    double int_tol_{1e-6};  // violation tolerance for integral feasibility check
    std::vector<std::unique_ptr<sep::Separator>> separators_;

    // Amortized separation: skip rounds to reduce overhead
    int32_t separation_interval_ = 1;  // 1 = every round (no skipping)
    int32_t threads_ = 1;              // the solve's thread budget
    int32_t max_cuts_per_sep_ = 3;     // max cuts per separator per round (0 = unlimited)
    int32_t max_cuts_per_round_ = 0;   // global cut cap per round (0 = unlimited)
    std::unordered_map<std::string, int32_t> per_separator_max_cuts_;
    bool submip_separation_ = true;  // SEC separation at sub-MIP root node
    double upper_bound_ = std::numeric_limits<double>::infinity();
    std::vector<double> warm_start_solution_;
    bool edge_elimination_enabled_ = true;
    bool edge_elimination_nodes_ = true;
    bool bounds_propagation_ = true;

    // LP-guided heuristic thresholds
    double heu_lpg_edge_threshold_ = 0.1;
    double heu_lpg_node_threshold_ = 0.5;
    double heu_lpg_lp_threshold_ = 0.1;
    double heu_lpg_seed_threshold_ = 0.3;

    // Labeling bounds for edge elimination and propagation
    std::vector<double> fwd_bounds_;               // f[v]: forward bounds (source → v)
    std::vector<double> bwd_bounds_;               // b[v]: backward bounds (v → target)
    double correction_ = 0.0;                      // profit(source) if tour, 0 if path
    std::vector<double> edge_lb_;                  // per-edge bounds, or empty (#144)
    std::vector<int8_t> fixed_y_;                  // optional fixed y-vars: -1 unset, 0/1 fixed
    std::vector<detail::Hyperplane> hyperplanes_;  // one z column each, after x and y

    // All-pairs labeling bounds, flat n×n: d(s,v) = (*all_pairs_)[s*n+v]. Shared,
    // not moved, between the separator callback (SPI) and the propagator
    // (Trigger B), which both read them for the whole solve (#156). Never null.
    std::shared_ptr<const std::vector<double>> all_pairs_ =
        std::make_shared<const std::vector<double>>();
    std::shared_ptr<std::atomic<bool>> interrupt_flag_;

    // RC fixing settings and statistics
    bool propagator_installed_ = false;
    bool obj_is_integer_ = false;  // true if all obj coefficients are integral
    int64_t labeling_max_queue_pops_ = 0;
    preprocess::LabelRows label_rows_ = preprocess::LabelRows::kFirst;
    preprocess::LabelingDeadline labeling_deadline_ = preprocess::kNoLabelingDeadline;
    RCFixingSettings rc_settings_;

    // Cut statistics, updated from the separator callback, which HiGHS's
    // parallel tree search calls from several workers at once: the counters
    // are atomic, and each round merges its per-separator stats and time
    // under cut_stats_mutex_.
    mutable std::mutex cut_stats_mutex_;
    mutable std::map<std::string, SeparatorStats> separator_stats_;
    mutable std::atomic<int64_t> total_cuts_ = 0;
    mutable std::atomic<int64_t> separation_rounds_ = 0;
    mutable std::atomic<int64_t> separation_calls_ = 0;  // total callback invocations
    // Wall-clock time of the separation phase: one span per round covering
    // support-graph/Gomory-Hu construction, the parallel separators, and cut
    // pooling (not the sum of per-family times, which would overcount the
    // parallelism). Summed over rounds, so over workers under parallel tree
    // search.
    mutable double separation_time_seconds_ = 0.0;

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

    // Callback statistics (single allocation shared with callback lambdas)
    std::shared_ptr<CallbackStats> stats_;

    // The propagator's state per tree worker, merged into stats_ at the end.
    std::shared_ptr<PropagatorStates> propagator_states_;
};

}  // namespace cptp
