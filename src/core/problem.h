#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "core/static_graph.h"

namespace cptp {

using Graph = static_graph;

struct Edge {
  int32_t tail;
  int32_t head;
};

/// Adjacency list entry: an edge index and the neighbor reachable via it.
struct AdjEntry {
  int32_t edge;
  int32_t neighbor;
};

/// CPTP instance: undirected graph with profits, costs, demands, capacity.
class Problem {
 public:
  Problem() = default;

  /// Build from raw data (copies from spans). edges[i] = {tail, head} with tail
  /// < head.
  void build(int32_t num_nodes, std::span<const Edge> edges,
             std::span<const double> edge_costs,
             std::span<const double> profits, std::span<const double> demands,
             double capacity, int32_t source = 0, int32_t target = 0);

  /// Build by moving owned vectors (avoids copies when caller no longer needs
  /// them).
  void build(int32_t num_nodes, std::span<const Edge> edges,
             std::vector<double> edge_costs, std::vector<double> profits,
             std::vector<double> demands, double capacity, int32_t source = 0,
             int32_t target = 0);

  // --- Accessors ---
  const Graph& graph() const { return graph_; }
  int32_t num_nodes() const { return num_nodes_; }
  int32_t num_edges() const { return num_edges_; }
  int32_t source() const { return source_; }
  int32_t target() const { return target_; }
  int32_t depot() const { return source_; }  // backward compat alias
  bool is_tour() const { return source_ == target_; }
  double capacity() const { return capacity_; }

  double edge_cost(int32_t e) const { return edge_costs_[e]; }
  double profit(int32_t v) const { return profits_[v]; }
  double demand(int32_t v) const { return demands_[v]; }

  const std::vector<double>& edge_costs() const { return edge_costs_; }
  const std::vector<double>& profits() const { return profits_; }
  const std::vector<double>& demands() const { return demands_; }

  std::string name;

 private:
  Graph graph_;
  int32_t num_nodes_ = 0;
  int32_t num_edges_ = 0;
  int32_t source_ = 0;
  int32_t target_ = 0;
  double capacity_ = 0.0;
  std::vector<double> edge_costs_;
  std::vector<double> profits_;
  std::vector<double> demands_;
};

}  // namespace cptp
