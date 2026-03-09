#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "core/digraph.h"
#include "core/dinitz.h"

namespace cptp {

/// Gomory-Hu tree via Gusfield's simplified algorithm.
/// Computes all pairwise min-cuts with n-1 max-flow calls.
class gomory_hu_tree {
 public:
  gomory_hu_tree(const digraph& g, const std::vector<double>& cap, int32_t root)
      : root_(root),
        n_(g.num_nodes()),
        parent_(n_, root),
        weight_(n_, 0.0),
        cut_side_(n_) {
    // Gusfield's algorithm: for each node i != root, compute
    // max-flow(i, parent[i]) and update tree edges.
    for (int32_t i = 0; i < n_; ++i) {
      if (i == root) continue;

      dinitz alg(g, cap, i, parent_[i]);
      alg.run();
      weight_[i] = alg.flow_value();

      // Store which nodes are on i's side (source side)
      cut_side_[i].resize(n_);
      for (int32_t v = 0; v < n_; ++v) {
        cut_side_[i][v] = alg.on_source_side(v);
      }

      // Update tree: for j > i where parent[j] == parent[i]
      // and j is on i's side, reparent j to i.
      for (int32_t j = i + 1; j < n_; ++j) {
        if (j == root) continue;
        if (parent_[j] == parent_[i] && cut_side_[i][j]) {
          parent_[j] = i;
        }
      }

      // Also: if weight[parent[i]] > weight[i] and parent[i]'s
      // cut puts i on source side, swap. (Gusfield's full update)
      // Already handled by the j > i reparenting above.
    }
  }

  struct CutInfo {
    double value;
    std::vector<bool> on_source_side;  // true = depot/root side
  };

  /// Lightweight reference to a stored cut (no copy of cut_side).
  struct CutRef {
    int32_t node;   // GH tree node whose edge defines this cut
    double weight;  // cut value
  };

  /// Return the min-cut separating target from root.
  /// Walks the tree path from target to root, finds the bottleneck.
  CutInfo min_cut(int32_t target) const {
    // Walk path from target to root, find minimum weight edge
    int32_t bottleneck = target;
    double min_weight = std::numeric_limits<double>::infinity();

    int32_t v = target;
    while (v != root_) {
      if (weight_[v] < min_weight) {
        min_weight = weight_[v];
        bottleneck = v;
      }
      v = parent_[v];
    }

    // The partition: nodes on bottleneck's side in its cut
    // cut_side_[bottleneck][v] = true means v was on source side
    // (bottleneck's side) of the max-flow bottleneck <-> parent[bottleneck].
    // We want on_source_side[v] = true for root side.
    // Since the max-flow was (bottleneck, parent[bottleneck]) with
    // bottleneck as source, on_source_side = !cut_side for root side.
    std::vector<bool> root_side(n_);
    for (int32_t u = 0; u < n_; ++u) {
      root_side[u] = !cut_side_[bottleneck][u];
    }

    return {min_weight, std::move(root_side)};
  }

  /// Return all cuts along the path from target to root.
  /// Each CutRef identifies a GH tree node whose stored cut_side
  /// defines the partition. Returns O(depth) candidates with zero
  /// additional max-flow cost.
  std::vector<CutRef> all_cuts_on_path(int32_t target) const {
    std::vector<CutRef> cuts;
    int32_t v = target;
    while (v != root_) {
      cuts.push_back({v, weight_[v]});
      v = parent_[v];
    }
    return cuts;
  }

  /// Get the cut partition for a given GH tree node.
  /// Returns true for nodes on the root/depot side (i.e., NOT on
  /// the node's source side in the original max-flow).
  bool on_root_side(int32_t cut_node, int32_t query_node) const {
    return !cut_side_[cut_node][query_node];
  }

  int32_t num_nodes() const { return n_; }

 private:
  int32_t root_;
  int32_t n_;
  std::vector<int32_t> parent_;
  std::vector<double> weight_;
  std::vector<std::vector<bool>> cut_side_;
};

}  // namespace cptp
