#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "core/problem.h"

namespace cptp::preprocess {

/// Label for the capacity-aware labeling.
/// Tracks net cost (edge costs - collected profits) and accumulated demand.
struct Label {
  double cost;    // net cost = sum(edge_costs) - sum(profits along path)
  double demand;  // accumulated demand along path
};

/// Capacity-aware labeling from a given root node, with custom edge costs and
/// node profits.
///
/// For each node v, computes the minimum net cost f[v] to reach v from root
/// via any (possibly non-elementary) path that respects vehicle capacity.
///
/// Net cost = sum of edge costs along path - sum of profits of visited nodes
/// (excluding root profit, which is always collected... see seed below).
///
/// Returns f[v] for each node. f[root] = -profits[root].
inline std::vector<double> labeling_from(const Problem& prob, int32_t root,
                                         std::span<const double> edge_costs,
                                         std::span<const double> profits,
                                         int64_t max_queue_pops = 0) {
  const int32_t n = prob.num_nodes();
  const auto& graph = prob.graph();
  const double Q = prob.capacity();
  constexpr double inf = std::numeric_limits<double>::infinity();

  std::vector<std::vector<Label>> labels(n);
  std::vector<double> best_cost(n, inf);

  double root_cost = -profits[root];
  double root_demand = prob.demand(root);
  labels[root].push_back({root_cost, root_demand});
  best_cost[root] = root_cost;

  struct QueueEntry {
    int32_t node;
    int32_t label_idx;
  };
  std::vector<QueueEntry> queue;
  queue.push_back({root, 0});

  size_t head = 0;
  while (head < queue.size()) {
    if (max_queue_pops > 0 && static_cast<int64_t>(head) >= max_queue_pops) {
      return {};  // budget exhausted — caller checks .empty()
    }
    auto [u, li] = queue[head++];

    if (li >= static_cast<int32_t>(labels[u].size())) continue;
    const Label& lbl = labels[u][li];

    for (auto e : graph.incident_edges(u)) {
      int32_t v = graph.other_endpoint(e, u);

      double new_demand = lbl.demand + prob.demand(v);
      if (new_demand > Q) continue;

      double new_cost = lbl.cost + edge_costs[e] - profits[v];

      auto& v_labels = labels[v];
      bool dominated = false;
      for (const auto& existing : v_labels) {
        if (existing.cost <= new_cost && existing.demand <= new_demand) {
          dominated = true;
          break;
        }
      }
      if (dominated) continue;

      std::erase_if(v_labels, [&](const Label& ex) {
        return new_cost <= ex.cost && new_demand <= ex.demand;
      });

      v_labels.push_back({new_cost, new_demand});

      if (new_cost < best_cost[v]) {
        best_cost[v] = new_cost;
      }

      queue.push_back({v, static_cast<int32_t>(v_labels.size()) - 1});
    }
  }

  return best_cost;
}

/// Capacity-aware labeling using the problem's own costs.
inline std::vector<double> labeling_from(const Problem& prob, int32_t root,
                                         int64_t max_queue_pops = 0) {
  return labeling_from(prob, root, prob.edge_costs(), prob.profits(),
                       max_queue_pops);
}

/// Forward labeling from a given source node (wrapper for labeling_from).
inline std::vector<double> forward_labeling(const Problem& prob,
                                            int32_t source,
                                            int64_t max_queue_pops = 0) {
  return labeling_from(prob, source, max_queue_pops);
}

/// Forward labeling from the problem's source node.
inline std::vector<double> forward_labeling(const Problem& prob,
                                            int64_t max_queue_pops = 0) {
  return labeling_from(prob, prob.source(), max_queue_pops);
}

/// Backward labeling from a given target node.
/// For undirected graphs, this is identical to forward_labeling(prob, target).
inline std::vector<double> backward_labeling(const Problem& prob,
                                             int32_t target,
                                             int64_t max_queue_pops = 0) {
  return labeling_from(prob, target, max_queue_pops);
}

/// Edge elimination using pre-computed forward/backward labeling bounds.
///
/// For each edge (u,v): lb = min(f[u]+c(u,v)+b[v], f[v]+c(u,v)+b[u]) +
/// correction If lb > upper_bound -> edge can be eliminated.
///
/// correction = profit(s) when s=t (tour, depot profit double-subtracted),
///              0 when s!=t (path).
///
/// Returns eliminated[e] = true for edges that can be fixed to 0.
inline std::vector<bool> edge_elimination(const Problem& prob,
                                          const std::vector<double>& f,
                                          const std::vector<double>& b,
                                          double upper_bound,
                                          double correction) {
  const auto& graph = prob.graph();
  const int32_t m = prob.num_edges();
  constexpr double inf = std::numeric_limits<double>::infinity();

  std::vector<bool> eliminated(m, false);

  if (upper_bound >= inf) return eliminated;

  for (auto e : graph.edges()) {
    int32_t u = graph.edge_source(e);
    int32_t v = graph.edge_target(e);

    if (f[u] >= inf || b[v] >= inf) {
      // Try both orientations for undirected
      if (f[v] >= inf || b[u] >= inf) {
        eliminated[e] = true;
        continue;
      }
    }

    // Lower bound on any tour/path using edge e:
    // min of both orientations (undirected edge)
    double lb1 = f[u] + prob.edge_cost(e) + b[v] + correction;
    double lb2 = f[v] + prob.edge_cost(e) + b[u] + correction;
    double lb = std::min(lb1, lb2);

    if (lb > upper_bound + 1e-6) {
      eliminated[e] = true;
    }
  }

  return eliminated;
}

/// Self-contained edge elimination: computes labeling bounds internally.
///
/// Tour: lb = f[u] + c(u,v) + f[v] + profit(depot)
/// Path: lb = f_s[u] + c(u,v) + f_t[v]
inline std::vector<bool> edge_elimination(const Problem& prob,
                                          double upper_bound) {
  const auto& graph = prob.graph();
  const int32_t m = prob.num_edges();
  constexpr double inf = std::numeric_limits<double>::infinity();

  std::vector<bool> eliminated(m, false);

  if (upper_bound >= inf) return eliminated;

  // Forward labeling from source
  auto f_s = labeling_from(prob, prob.source());

  if (prob.is_tour()) {
    // Tour: f_s[u] + c(u,v) + f_s[v] + profit(depot)
    // (undirected: backward = forward, correct for double depot subtraction)
    double depot_profit = prob.profit(prob.source());

    for (auto e : graph.edges()) {
      int32_t u = graph.edge_source(e);
      int32_t v = graph.edge_target(e);

      if (f_s[u] >= inf || f_s[v] >= inf) {
        eliminated[e] = true;
        continue;
      }

      double lb = f_s[u] + prob.edge_cost(e) + f_s[v] + depot_profit;
      if (lb > upper_bound + 1e-6) {
        eliminated[e] = true;
      }
    }
  } else {
    // Path: need labeling from both source and target
    auto f_t = labeling_from(prob, prob.target());

    for (auto e : graph.edges()) {
      int32_t u = graph.edge_source(e);
      int32_t v = graph.edge_target(e);

      // Try both orientations: source->...->u->v->...->target and reverse
      double lb1 = (f_s[u] < inf && f_t[v] < inf)
                       ? f_s[u] + prob.edge_cost(e) + f_t[v]
                       : inf;
      double lb2 = (f_s[v] < inf && f_t[u] < inf)
                       ? f_s[v] + prob.edge_cost(e) + f_t[u]
                       : inf;
      double lb = std::min(lb1, lb2);

      if (lb >= inf || lb > upper_bound + 1e-6) {
        eliminated[e] = true;
      }
    }
  }

  return eliminated;
}

}  // namespace cptp::preprocess
