#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace cptp {

/// Per-separator statistics collected during branch-and-cut.
struct SeparatorStats {
  int64_t cuts_added = 0;     // total cuts added across all rounds
  int64_t rounds_called = 0;  // number of separation rounds
  double time_seconds = 0.0;  // cumulative wall-clock time in separator
};

struct SolveResult {
  enum class Status {
    Optimal,
    Feasible,
    Infeasible,
    Unbounded,
    TimeLimit,
    Error
  };

  Status status = Status::Error;
  double objective = 0.0;
  double bound = 0.0;
  double gap = 1.0;
  double time_seconds = 0.0;
  int64_t nodes = 0;
  int64_t simplex_iterations = -1;

  std::vector<int32_t> tour;       // ordered node indices
  std::vector<int32_t> tour_arcs;  // arc indices in tour

  /// Cut statistics per separator (keyed by separator name).
  std::map<std::string, SeparatorStats> separator_stats;
  int64_t total_cuts = 0;
  int64_t separation_rounds = 0;

  /// Per-phase wall-clock times (seconds), accumulated on the callback thread.
  /// Callbacks run single-threaded; their internal parallelism (parallel
  /// separators, parallel labeling) is measured as one span per call, so these
  /// credit the parallelism rather than summing concurrent work.
  double separation_time_seconds = 0.0;
  double propagator_time_seconds = 0.0;
  double rc_time_seconds = 0.0;  // subset of propagator_time_seconds (runs inside it)
  double heuristic_time_seconds = 0.0;  // LP-guided B&C callback only (not warm-start)

  /// Domain-propagation fixings (edges from sweep/chain, plus node fixings).
  int64_t sweep_fixings = 0;
  int64_t chain_fixings = 0;
  int64_t sweep_node_fixings = 0;
  int64_t chain_node_fixings = 0;

  /// Reduced-cost fixings: variables fixed to 0 and to 1.
  int64_t rc_fix0_count = 0;
  int64_t rc_fix1_count = 0;

  bool is_optimal() const { return status == Status::Optimal; }
  bool has_solution() const {
    return status == Status::Optimal || status == Status::Feasible ||
           status == Status::TimeLimit;
  }
};

}  // namespace cptp
