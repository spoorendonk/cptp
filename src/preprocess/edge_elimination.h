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

/// Capacity-aware 2-cycle elimination labeling from the depot.
///
/// For each node v, computes the minimum net cost f[v] to reach v from depot
/// via any path that respects: (1) vehicle capacity, (2) 2-cycle elimination
/// (no immediate return to predecessor).
///
/// Net cost = sum of edge costs along path - sum of profits of visited nodes
/// (excluding depot profit, which is always collected).
///
/// Returns f[v] for each node. f[depot] = -profit(depot).
inline std::vector<double> forward_labeling(const Problem& prob) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();
    constexpr double inf = std::numeric_limits<double>::infinity();

    // For each node, maintain a set of non-dominated labels.
    // Label dominance: (cost1, demand1) dominates (cost2, demand2) iff
    //   cost1 <= cost2 AND demand1 <= demand2.
    std::vector<std::vector<Label>> labels(n);

    // Best (minimum) cost to reach each node (for quick lower bound queries).
    std::vector<double> best_cost(n, inf);

    // Seed: depot with cost = -profit(depot), demand = demand(depot), prev = -1
    double depot_cost = -prob.profit(depot);
    double depot_demand = prob.demand(depot);
    labels[depot].push_back({depot_cost, depot_demand, -1});
    best_cost[depot] = depot_cost;

    // BFS/label-setting with a queue of (node, label_index) to propagate.
    // Since costs can be negative (from profits), we cannot use Dijkstra.
    // Use FIFO label-correcting approach with dominance checks.
    //
    // To bound complexity, we limit labels per node. For edge elimination
    // purposes, we only need good lower bounds, not exact shortest paths.
    constexpr int32_t kMaxLabelsPerNode = 50;

    struct QueueEntry {
        int32_t node;
        int32_t label_idx;
    };
    std::vector<QueueEntry> queue;
    queue.push_back({depot, 0});

    size_t head = 0;
    while (head < queue.size()) {
        auto [u, li] = queue[head++];

        // Check if this label still exists (may have been dominated)
        if (li >= static_cast<int32_t>(labels[u].size())) continue;
        const Label& lbl = labels[u][li];

        // Propagate along all incident edges
        for (auto e : graph.incident_edges(u)) {
            int32_t v = graph.other_endpoint(e, u);

            // 2-cycle elimination: don't go back to predecessor
            if (v == lbl.prev) continue;

            double new_demand = lbl.demand + prob.demand(v);
            if (new_demand > Q) continue;  // capacity violated

            double new_cost = lbl.cost + prob.edge_cost(e) - prob.profit(v);

            // Quick prune: if new_cost is worse than best for v, skip
            // (not exact dominance, but avoids many labels)
            // We still add if demand is significantly lower than existing.

            // Dominance check against existing labels at v
            auto& v_labels = labels[v];
            bool dominated = false;
            for (const auto& existing : v_labels) {
                if (existing.cost <= new_cost && existing.demand <= new_demand) {
                    dominated = true;
                    break;
                }
            }
            if (dominated) continue;

            // Remove labels dominated by the new one
            std::erase_if(v_labels, [&](const Label& ex) {
                return new_cost <= ex.cost && new_demand <= ex.demand;
            });

            if (static_cast<int32_t>(v_labels.size()) >= kMaxLabelsPerNode) {
                // Only add if better cost than worst existing label
                auto worst = std::max_element(v_labels.begin(), v_labels.end(),
                    [](const Label& a, const Label& b) { return a.cost < b.cost; });
                if (new_cost >= worst->cost) continue;
                // Replace worst
                *worst = {new_cost, new_demand, u};
            } else {
                v_labels.push_back({new_cost, new_demand, u});
            }

            int32_t new_idx = static_cast<int32_t>(v_labels.size()) - 1;
            // Find the actual index (might have replaced worst)
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

/// Edge elimination using capacity-aware 2-cycle labeling bounds.
///
/// For each edge (u,v) with cost c:
///   lb(u,v) = f[u] + c - profit(v) + f[v]  (use edge to reach v from u's side)
///   but also = f[v] + c - profit(u) + f[u]  (symmetric, use edge to reach u from v's side)
///   Tighter: lb = min(f[u] + c + g[v], f[v] + c + g[u])
///   where g[v] = min cost to reach depot FROM v (= f[v] for undirected graphs).
///
/// Actually, for undirected the lower bound on a round-trip using edge (u,v) is:
///   lb = f[u] + c + f[v]
///   where f includes the -profit terms already (net cost to reach node from depot).
///   But we need to add back the depot profit since f[depot] = -profit(depot) and
///   the objective includes depot cost only once.
///
/// If lb > UB → edge can be eliminated.
///
/// Also eliminates nodes: if all incident edges are eliminated, fix y_v = 0.
///
/// Returns eliminated[e] = true for edges that can be fixed to 0.
inline std::vector<bool> edge_elimination(const Problem& prob, double upper_bound) {
    const auto& graph = prob.graph();
    const int32_t m = prob.num_edges();
    constexpr double inf = std::numeric_limits<double>::infinity();

    std::vector<bool> eliminated(m, false);

    if (upper_bound >= inf) return eliminated;  // no bound, can't eliminate

    // Forward labeling: f[v] = min net cost to reach v from depot
    auto f = forward_labeling(prob);

    // For undirected graphs, backward = forward, so g[v] = f[v].
    // But the round-trip bound needs care:
    //   A tour using edge (u,v) costs at least:
    //     f[u] + edge_cost(u,v) + f[v]
    //   But this double-counts depot profit (included in both f[u] and f[v] paths
    //   only if they share nodes). Actually, f[v] is the net cost depot→...→v,
    //   so f[u] + c(u,v) + f[v] estimates depot→...→u→v→...→depot.
    //   The depot cost -profit(depot) is only in the initial seed, and the two
    //   paths share only the depot node. So the round-trip cost is:
    //     f[u] + edge_cost(u,v) + f[v] + profit(depot)
    //   because f includes -profit(depot) in the initial seed for BOTH directions,
    //   but depot profit should only be counted once in the objective.
    //
    //   Wait — f[v] represents "cost to reach v starting from depot". The return
    //   path f[v] represents "cost from depot to v" = "cost from v to depot" (undirected).
    //   So a round trip using edge (u,v):
    //     depot →...→ u → v →...→ depot
    //   costs: f[u] + c(u,v) - profit(v) + (cost from v to depot without profit(v))
    //   But f[v] already includes -profit(v)... this is getting circular.
    //
    //   Simpler: f[v] = min net cost path depot→v = sum(edge costs) - sum(node profits except what's collected).
    //   The labeling subtracts profit(v) when arriving at v.
    //   So f[u] = net cost depot→...→u (profits of all nodes on path subtracted).
    //   Similarly f[v] = net cost depot→...→v.
    //   Round trip through (u,v): depot→...→u→v→...→depot
    //   Lower bound = f[u] + edge_cost(u,v) - profit(v) + f_return[v]
    //   But f_return[v] = f[v] + profit(depot) because the return path from v to
    //   depot also starts with depot seed = -profit(depot), and depot profit is
    //   collected in both directions.
    //
    //   Let's just use: lb = f[u] + edge_cost(u,v) + f[v] + profit(depot)
    //   This is correct because both f[u] and f[v] subtract profit(depot) in
    //   their seed, but depot appears once in the tour.

    double depot_profit = prob.profit(prob.depot());

    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);

        if (f[u] >= inf || f[v] >= inf) {
            eliminated[e] = true;
            continue;
        }

        // Lower bound on any tour using edge e:
        // depot→...→u→v→...→depot
        // = f[u] + c(u,v) + f[v] + profit(depot)
        // (adding back depot profit because it's double-subtracted in f[u] and f[v])
        double lb = f[u] + prob.edge_cost(e) + f[v] + depot_profit;

        if (lb > upper_bound + 1e-6) {
            eliminated[e] = true;
        }
    }

    return eliminated;
}

}  // namespace cptp::preprocess
