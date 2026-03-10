#pragma once

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace cptp {

/// Simple undirected graph with CSR storage.
class static_graph {
 public:
  struct Edge {
    int32_t u, v;  // u < v
  };

  static_graph() = default;

  static_graph(int32_t n, std::span<const Edge> edges)
      : num_nodes_(n), num_edges_(static_cast<int32_t>(edges.size())) {
    endpoints_.resize(2 * num_edges_);
    for (int32_t i = 0; i < num_edges_; ++i) {
      if (edges[i].u >= edges[i].v)
        throw std::invalid_argument("static_graph: edge u must be less than v");
      endpoints_[2 * i] = edges[i].u;
      endpoints_[2 * i + 1] = edges[i].v;
    }

    // Build incidence lists (CSR)
    // Count incident edges per node
    std::vector<int32_t> degree(n, 0);
    for (int32_t e = 0; e < num_edges_; ++e) {
      degree[endpoints_[2 * e]]++;
      degree[endpoints_[2 * e + 1]]++;
    }

    incidence_begin_.resize(n + 1);
    incidence_begin_[0] = 0;
    for (int32_t i = 0; i < n; ++i) {
      incidence_begin_[i + 1] = incidence_begin_[i] + degree[i];
    }

    incidence_.resize(incidence_begin_[n]);
    std::vector<int32_t> pos(n, 0);
    for (int32_t e = 0; e < num_edges_; ++e) {
      int32_t u = endpoints_[2 * e];
      int32_t v = endpoints_[2 * e + 1];
      incidence_[incidence_begin_[u] + pos[u]++] = e;
      incidence_[incidence_begin_[v] + pos[v]++] = e;
    }
  }

  int32_t num_nodes() const { return num_nodes_; }
  int32_t num_edges() const { return num_edges_; }

  /// Iterate all edge indices [0, num_edges).
  struct edge_range {
    int32_t n;
    struct iterator {
      int32_t i;
      int32_t operator*() const { return i; }
      iterator& operator++() {
        ++i;
        return *this;
      }
      bool operator!=(const iterator& o) const { return i != o.i; }
    };
    iterator begin() const { return {0}; }
    iterator end() const { return {n}; }
  };
  edge_range edges() const { return {num_edges_}; }

  /// Edges incident to node v.
  std::span<const int32_t> incident_edges(int32_t v) const {
    return {incidence_.data() + incidence_begin_[v],
            incidence_.data() + incidence_begin_[v + 1]};
  }

  int32_t edge_source(int32_t e) const { return endpoints_[2 * e]; }
  int32_t edge_target(int32_t e) const { return endpoints_[2 * e + 1]; }

  int32_t other_endpoint(int32_t e, int32_t v) const {
    int32_t u = endpoints_[2 * e];
    int32_t w = endpoints_[2 * e + 1];
    return (v == u) ? w : u;
  }

 private:
  int32_t num_nodes_ = 0;
  int32_t num_edges_ = 0;
  std::vector<int32_t> endpoints_;        // [2*e, 2*e+1] = {source, target}
  std::vector<int32_t> incidence_begin_;  // CSR offsets, size n+1
  std::vector<int32_t> incidence_;        // edge indices
};

}  // namespace cptp
