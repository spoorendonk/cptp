#pragma once

#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "core/problem.h"

namespace cptp::preprocess {

/// Solver-independent domain propagation based on labeling bounds.
///
/// Given pre-computed forward/backward labeling bounds and the current
/// upper bound, determines which edge/node variables can be fixed to 0.
///
/// Two propagation mechanisms:
///   - Trigger A (sweep): when UB improves, scan all edges.
///   - Trigger B (chain): when an edge is fixed to 1, infer fixings for others.
///
/// Usage:
///   BoundPropagator prop(prob, fwd, bwd, correction);
///   auto edges = prop.sweep(upper_bound, col_upper);
///   auto nodes = prop.sweep_nodes(col_upper, y_offset);
///   auto chain = prop.propagate_fixed_edge(e, upper_bound, col_upper);
class BoundPropagator {
 public:
  BoundPropagator(const Problem& prob, std::vector<double> fwd_bounds,
                  std::vector<double> bwd_bounds, double correction);

  /// Set all-pairs labeling bounds for stronger Trigger B propagation.
  /// dist is a flat n*n matrix: d(s,v) = dist[s*n + v].
  void set_all_pairs_bounds(std::vector<double> dist);

  bool has_all_pairs_bounds() const { return !all_pairs_.empty(); }

  /// Trigger A: full sweep — return edge indices that can be fixed to 0.
  /// col_upper[e] is the current upper bound of edge variable e.
  /// Only edges with col_upper[e] > 0.5 (i.e., not already fixed) are checked.
  std::vector<int32_t> sweep(double upper_bound,
                             std::span<const double> col_upper) const;

  /// After a sweep, find nodes whose all incident edges are fixed to 0.
  /// col_upper covers edges [0, m) and nodes [m, m+n).
  /// y_offset is the column offset for node variables (typically = num_edges).
  /// Returns node variable indices (y_offset + i) that can be fixed to 0.
  std::vector<int32_t> sweep_nodes(std::span<const double> col_upper,
                                   int32_t y_offset) const;

  /// Trigger B: given edge e fixed to 1, return other edges that can be
  /// fixed to 0.  col_upper[e] is current upper bound of edge variable e.
  std::vector<int32_t> propagate_fixed_edge(
      int32_t edge, double upper_bound,
      std::span<const double> col_upper) const;

  const std::vector<double>& fwd_bounds() const { return fwd_bounds_; }
  const std::vector<double>& bwd_bounds() const { return bwd_bounds_; }
  double correction() const { return correction_; }

 private:
  const Problem& prob_;
  std::vector<double> fwd_bounds_;
  std::vector<double> bwd_bounds_;
  double correction_;
  std::vector<double> all_pairs_;
  std::vector<double> edge_costs_;
  std::vector<double> profits_;
  std::vector<std::vector<AdjEntry>> adjacency_;
};

}  // namespace cptp::preprocess
