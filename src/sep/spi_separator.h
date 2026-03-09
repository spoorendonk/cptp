#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Shortest Path Inequality (SPI) separator.
///
/// Uses all-pairs shortest path bounds to derive infeasibility cuts.
/// Given an upper bound UB and shortest-path distances d(u,v), if the
/// cheapest tour/path visiting a set of nodes S exceeds UB, the nodes
/// in S cannot all be selected simultaneously.
///
/// Pair:    y_i + y_j <= 1           (any tour visiting both costs > UB)
/// Triplet: y_i + y_j + y_k <= 2    (any tour visiting all three costs > UB)
///
/// The lower bound for a pair {i,j} on a tour is:
///   lb = min over orderings of d(depot,i) + d(i,j) + d(j,depot)
///        + profit(depot) + profit(i) + profit(j)
/// The junction-node profit correction accounts for profits being
/// double-subtracted when concatenating shortest-path segments.
class SPISeparator : public Separator {
 public:
  std::string name() const override { return "SPI"; }
  std::vector<Cut> separate(const SeparationContext& ctx) override;

  void set_max_dp_size(int32_t v) { max_dp_size_ = v; }
  void set_max_seeds(int32_t v) { max_seeds_ = v; }
  void set_max_grow_candidates(int32_t v) { max_grow_candidates_ = v; }

  int32_t max_dp_size() const { return max_dp_size_; }
  int32_t max_seeds() const { return max_seeds_; }
  int32_t max_grow_candidates() const { return max_grow_candidates_; }

 private:
  int32_t max_dp_size_ = 15;
  int32_t max_seeds_ = 50;
  int32_t max_grow_candidates_ = 40;
};

}  // namespace cptp::sep
