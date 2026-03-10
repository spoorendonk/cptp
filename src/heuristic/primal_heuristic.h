#pragma once

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <vector>

#include "core/problem.h"
#include "sep/separation_oracle.h"

namespace cptp::heuristic {

/// Subset of the original problem's edges/nodes, built from LP relaxation
/// values.
struct ReducedGraph {
  std::vector<bool>
      edge_active;  // size m: which edges are in the reduced graph
  std::vector<bool> node_active;  // size n: which nodes are candidates
};

inline constexpr int32_t kDefaultLsMaxIterPerStart = 1000;

namespace detail {

inline constexpr double kNoSolution = std::numeric_limits<double>::infinity();
inline constexpr double kLsEpsImprove = 1e-9;
inline constexpr double kLsEpsTie = 1e-12;
inline constexpr double kBigEdgeCost = 1e15;
inline constexpr int32_t kIlsRounds = 8;
inline constexpr int32_t kKickTrials = 16;

struct TourCandidate {
  std::vector<int32_t> tour;
  double objective = kNoSolution;
};

/// Find edge index between u and v, or -1 if not found.
/// When edge_active is non-empty, inactive edges return -1.
inline int32_t find_edge(const Graph& g, int32_t u, int32_t v,
                         const std::vector<bool>& edge_active = {}) {
  for (int32_t e : g.incident_edges(u)) {
    if (!edge_active.empty() && !edge_active[e]) continue;
    if (g.other_endpoint(e, u) == v) return e;
  }
  return -1;
}

inline double edge_cost(const Problem& prob, int32_t u, int32_t v,
                        const std::vector<bool>& edge_active = {}) {
  int32_t e = find_edge(prob.graph(), u, v, edge_active);
  return (e >= 0) ? prob.edge_cost(e) : std::numeric_limits<double>::max();
}

/// Compute tour/path objective: sum(edge costs) - sum(profits).
double tour_objective(const Problem& prob, const std::vector<int32_t>& tour);

/// Greedy insertion: insert customers in given order at their cheapest
/// position. When edge_active is non-empty, only active edges are considered.
void greedy_insert(const Problem& prob, std::vector<int32_t>& tour,
                   std::vector<bool>& in_tour, double& remaining_cap,
                   const std::vector<int32_t>& order,
                   const std::vector<bool>& edge_active = {});

/// Run local search neighborhoods until no improvement or iteration limit.
/// When edge_active is non-empty, only active edges are considered.
int local_search(const Problem& prob, std::vector<int32_t>& tour,
                 std::vector<bool>& in_tour, double& remaining_cap,
                 const std::vector<bool>& edge_active = {},
                 int max_iter_per_start = kDefaultLsMaxIterPerStart);

inline bool seed_uses_active_edges(const Problem& prob,
                                   std::span<const int32_t> tour,
                                   const std::vector<bool>& edge_active = {}) {
  for (size_t i = 0; i + 1 < tour.size(); ++i) {
    if (edge_cost(prob, tour[i], tour[i + 1], edge_active) > kBigEdgeCost) {
      return false;
    }
  }
  return true;
}

bool init_seed_state(const Problem& prob, std::span<const int32_t> tour,
                     std::vector<bool>& in_tour, double& remaining_cap);

/// SplitMix64 for deterministic index mixing (no RNG state, no time seeding).
inline uint64_t mix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}

inline uint64_t deterministic_kick_key(int32_t start_index, int32_t round_index,
                                       int32_t trial, uint64_t salt) {
  uint64_t key = static_cast<uint64_t>(start_index + 1);
  key ^= static_cast<uint64_t>(round_index + 1) * 0x9e3779b97f4a7c15ull;
  key ^= static_cast<uint64_t>(trial + 1) * 0xbf58476d1ce4e5b9ull;
  key ^= salt;
  return mix64(key);
}

bool perturb_swap_kick(const Problem& prob, std::vector<int32_t>& tour,
                       std::vector<bool>& in_tour, double& remaining_cap,
                       const std::vector<bool>& edge_active,
                       int32_t start_index, int32_t round_index);

bool perturb_relocate_kick(const Problem& prob, std::vector<int32_t>& tour,
                           const std::vector<bool>& edge_active,
                           int32_t start_index, int32_t round_index);

inline bool perturb_deterministic(const Problem& prob,
                                  std::vector<int32_t>& tour,
                                  std::vector<bool>& in_tour,
                                  double& remaining_cap,
                                  const std::vector<bool>& edge_active,
                                  int32_t start_index, int32_t round_index) {
  if (perturb_swap_kick(prob, tour, in_tour, remaining_cap, edge_active,
                        start_index, round_index)) {
    return true;
  }
  return perturb_relocate_kick(prob, tour, edge_active, start_index,
                               round_index);
}

/// Single restart: construct + local search, return (tour, objective).
std::pair<std::vector<int32_t>, double> single_restart(
    const Problem& prob, const std::vector<int32_t>& customers,
    const std::vector<int32_t>& order,
    const std::vector<bool>& edge_active = {});

/// Deterministic fallback constructor used when restart search finds no
/// feasible tour/path.
std::vector<int32_t> fallback_feasible(
    const Problem& prob, const std::vector<bool>& edge_active = {});

/// Enumerate all 1-customer / 2-customer seeds and return sorted candidates.
std::vector<TourCandidate> enumerate_small_seed_candidates(
    const Problem& prob, const std::vector<int32_t>& customers,
    const std::vector<bool>& edge_active = {});

/// Convert a tour (node sequence) to a MIP solution vector.
void validate_tour_or_throw(const Problem& prob,
                            std::span<const int32_t> tour);

std::vector<double> tour_to_solution(const Problem& prob,
                                     const std::vector<int32_t>& tour);

}  // namespace detail

struct HeuristicResult {
  std::vector<double> col_values;  // size (num_edges + num_nodes)
  double objective;                // minimization objective (cost - profit)
};

struct WarmStartProgressSnapshot {
  int32_t starts_total = 0;
  int32_t starts_done = 0;
  int64_t ls_iterations_total = 0;
  int64_t ub_improvements = 0;
  double best_objective = detail::kNoSolution;
  bool final = false;
};

struct WarmStartProgressOptions {
  int64_t report_every_starts = 8;
  bool report_on_ub_improvement = true;
};

using WarmStartProgressCallback =
    std::function<void(const WarmStartProgressSnapshot&)>;

struct ConstructionPool {
  std::vector<detail::TourCandidate> candidates;  // sorted by objective
  std::vector<bool> edge_active;                  // optional reduced graph
};

// --- Graph reduction strategies ---

/// Strategy A: LP-value threshold.
ReducedGraph reduce_lp_threshold(const Problem& prob,
                                 std::span<const double> x_lp,
                                 std::span<const double> y_lp,
                                 double edge_threshold = 0.1,
                                 double node_threshold = 0.5);

/// Strategy B: RINS-style (when incumbent exists).
ReducedGraph reduce_rins(const Problem& prob, std::span<const double> x_lp,
                         std::span<const double> y_lp,
                         std::span<const double> incumbent,
                         double lp_threshold = 0.1);

/// Strategy C: Neighborhood expansion.
ReducedGraph reduce_neighborhood(const Problem& prob,
                                 std::span<const double> x_lp,
                                 std::span<const double> y_lp,
                                 double seed_threshold = 0.3);

// --- Public heuristic functions ---

std::vector<bool> build_edge_activity_from_bounds(
    const Problem& prob, std::span<const double> fwd_bounds,
    std::span<const double> bwd_bounds, double correction, double upper_bound);

ConstructionPool build_construction_pool(
    const Problem& prob, int max_candidates = 0,
    std::span<const double> fwd_bounds = {},
    std::span<const double> bwd_bounds = {}, double correction = 0.0,
    double upper_bound = std::numeric_limits<double>::infinity());

inline HeuristicResult best_construction_solution(
    const Problem& prob, const ConstructionPool& pool) {
  if (pool.candidates.empty()) {
    auto fallback = detail::fallback_feasible(prob, pool.edge_active);
    if (fallback.empty()) return {{}, detail::kNoSolution};
    const double obj = detail::tour_objective(prob, fallback);
    return {detail::tour_to_solution(prob, fallback), obj};
  }

  const auto& best = pool.candidates.front();
  return {detail::tour_to_solution(prob, best.tour), best.objective};
}

HeuristicResult run_local_search_from_pool(
    const Problem& prob, const ConstructionPool& pool, int num_starts,
    int ls_max_iter_per_start = kDefaultLsMaxIterPerStart,
    int32_t max_workers = 0, const std::vector<bool>& edge_active_override = {},
    const WarmStartProgressOptions* progress_opts = nullptr,
    const WarmStartProgressCallback& progress_cb = {},
    bool reuse_candidates_for_starts = false);

/// Build an initial solution via deterministic seed construction + parallel
/// local search over top seeds.
HeuristicResult build_initial_solution(
    const Problem& prob, int num_restarts = 50, double time_budget_ms = 0.0,
    std::span<const double> fwd_bounds = {},
    std::span<const double> bwd_bounds = {}, double correction = 0.0,
    double upper_bound = std::numeric_limits<double>::infinity(),
    int32_t max_workers = 0);

/// Compatibility wrapper: build_warm_start delegates to build_initial_solution.
inline HeuristicResult build_warm_start(
    const Problem& prob, int num_restarts = 50, double time_budget_ms = 0.0,
    std::span<const double> fwd_bounds = {},
    std::span<const double> bwd_bounds = {}, double correction = 0.0,
    double upper_bound = std::numeric_limits<double>::infinity(),
    int32_t max_workers = 0) {
  return build_initial_solution(prob, num_restarts, time_budget_ms, fwd_bounds,
                                bwd_bounds, correction, upper_bound,
                                max_workers);
}

/// LP-guided primal heuristic: builds reduced graphs from LP relaxation
/// values and runs combinatorial construction + local search on them.
HeuristicResult lp_guided_heuristic(
    const Problem& prob, std::span<const double> x_lp,
    std::span<const double> y_lp, std::span<const double> incumbent,
    double time_budget_ms = 20.0, int strategy = 0, int max_restarts = 0,
    uint32_t restart_seed = 0, double lpg_edge_threshold = 0.1,
    double lpg_node_threshold = 0.5, double lpg_lp_threshold = 0.1,
    double lpg_seed_threshold = 0.3);

}  // namespace cptp::heuristic
