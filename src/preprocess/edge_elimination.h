#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/problem.h"

namespace cptp::preprocess {

/// Label for the forward/backward ESPPRC-style labeling.
/// Tracks net cost (edge costs - collected profits) and accumulated demand.
struct Label {
    double cost;     // net cost = sum(edge_costs) - sum(profits along path)
    double demand;   // accumulated demand along path
    int32_t prev;    // previous node (for 2-cycle elimination)
};

/// Capacity-aware 2-cycle elimination labeling from a given root node.
///
/// For each node v, computes the minimum net cost f[v] to reach v from root
/// via any path that respects: (1) vehicle capacity, (2) 2-cycle elimination
/// (no immediate return to predecessor).
///
/// Net cost = sum of edge costs along path - sum of profits of visited nodes
/// (excluding root profit, which is always collected... see seed below).
///
/// Returns f[v] for each node. f[root] = -profit(root).
inline std::vector<double> labeling_from(const Problem& prob, int32_t root) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    const double Q = prob.capacity();
    constexpr double inf = std::numeric_limits<double>::infinity();

    // For each node, maintain a set of non-dominated labels.
    // Label dominance: (cost1, demand1) dominates (cost2, demand2) iff
    //   cost1 <= cost2 AND demand1 <= demand2.
    std::vector<std::vector<Label>> labels(n);

    // Best (minimum) cost to reach each node (for quick lower bound queries).
    std::vector<double> best_cost(n, inf);

    // Seed: root with cost = -profit(root), demand = demand(root), prev = -1
    double root_cost = -prob.profit(root);
    double root_demand = prob.demand(root);
    labels[root].push_back({root_cost, root_demand, -1});
    best_cost[root] = root_cost;

    // BFS/label-setting with a queue of (node, label_index) to propagate.
    // Since costs can be negative (from profits), we cannot use Dijkstra.
    // Use FIFO label-correcting approach with dominance checks.
    constexpr int32_t kMaxLabelsPerNode = 50;

    struct QueueEntry {
        int32_t node;
        int32_t label_idx;
    };
    std::vector<QueueEntry> queue;
    queue.push_back({root, 0});

    size_t head = 0;
    while (head < queue.size()) {
        auto [u, li] = queue[head++];

        if (li >= static_cast<int32_t>(labels[u].size())) continue;
        const Label& lbl = labels[u][li];

        for (auto e : graph.incident_edges(u)) {
            int32_t v = graph.other_endpoint(e, u);

            if (v == lbl.prev) continue;

            double new_demand = lbl.demand + prob.demand(v);
            if (new_demand > Q) continue;

            double new_cost = lbl.cost + prob.edge_cost(e) - prob.profit(v);

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

            if (static_cast<int32_t>(v_labels.size()) >= kMaxLabelsPerNode) {
                auto worst = std::max_element(v_labels.begin(), v_labels.end(),
                    [](const Label& a, const Label& b) { return a.cost < b.cost; });
                if (new_cost >= worst->cost) continue;
                *worst = {new_cost, new_demand, u};
            } else {
                v_labels.push_back({new_cost, new_demand, u});
            }

            int32_t new_idx = static_cast<int32_t>(v_labels.size()) - 1;
            for (int32_t i = 0; i < static_cast<int32_t>(v_labels.size()); ++i) {
                if (v_labels[i].cost == new_cost && v_labels[i].demand == new_demand
                    && v_labels[i].prev == u) {
                    new_idx = i;
                    break;
                }
            }

            if (new_cost < best_cost[v]) {
                best_cost[v] = new_cost;
            }

            queue.push_back({v, new_idx});
        }
    }

    return best_cost;
}

/// Backward-compatible wrapper: forward labeling from depot.
inline std::vector<double> forward_labeling(const Problem& prob) {
    return labeling_from(prob, prob.source());
}

/// Edge elimination using capacity-aware 2-cycle labeling bounds.
///
/// Tour: lb = f[u] + c(u,v) + f[v] + profit(depot)
///   (profit(depot) corrects double-subtraction of depot profit)
///
/// Path: lb = f_s[u] + c(u,v) + f_t[v]
///   (no double-counting since source != target, but both f_s and f_t
///    subtract their root profits. We need to add back one overlap correction
///    only if source == target, which is the tour case.)
///
/// If lb > UB -> edge can be eliminated.
inline std::vector<bool> edge_elimination(const Problem& prob, double upper_bound) {
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

            // Try both orientations: source→...→u→v→...→target and source→...→v→u→...→target
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
