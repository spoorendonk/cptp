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

  bool is_optimal() const { return status == Status::Optimal; }
  bool has_solution() const {
    return status == Status::Optimal || status == Status::Feasible ||
           status == Status::TimeLimit;
  }
};

}  // namespace cptp
