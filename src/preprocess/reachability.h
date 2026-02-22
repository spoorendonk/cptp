#pragma once

#include <limits>
#include <queue>
#include <vector>

#include "core/problem.h"

namespace rcspp::preprocess {

/// Dijkstra from a given root with weight = demand(target_node).
/// Returns dist[v] = minimum cumulative demand on any path from root to v.
inline std::vector<double> demand_dijkstra(const Problem& prob, int32_t root) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    constexpr double inf = std::numeric_limits<double>::infinity();

    std::vector<double> dist(n, inf);
    dist[root] = 0.0;

    using Entry = std::pair<double, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
    pq.emplace(0.0, root);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;

        for (auto e : graph.incident_edges(u)) {
            int32_t v = graph.other_endpoint(e, u);
            double nd = d + prob.demand(v);
            if (nd < dist[v]) {
                dist[v] = nd;
                pq.emplace(nd, v);
            }
        }
    }

    return dist;
}

/// Demand-reachability preprocessing.
/// Tour (source == target): round-trip check 2*dist[v] - demand(v) <= Q.
/// Path (source != target): one-way check dist_s[v] + dist_t[v] - demand(v) <= Q.
inline std::vector<bool> demand_reachability(const Problem& prob) {
    const int32_t n = prob.num_nodes();
    const double Q = prob.capacity();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    constexpr double inf = std::numeric_limits<double>::infinity();

    auto dist_s = demand_dijkstra(prob, source);

    std::vector<bool> reachable(n, false);

    if (prob.is_tour()) {
        // Tour: round-trip check
        for (int32_t v = 0; v < n; ++v) {
            if (dist_s[v] < inf) {
                reachable[v] = (2.0 * dist_s[v] - prob.demand(v) <= Q);
            }
        }
    } else {
        // Path: Dijkstra from both source and target
        auto dist_t = demand_dijkstra(prob, target);
        for (int32_t v = 0; v < n; ++v) {
            if (dist_s[v] < inf && dist_t[v] < inf) {
                reachable[v] = (dist_s[v] + dist_t[v] - prob.demand(v) <= Q);
            }
        }
    }

    // Source and target are always reachable
    reachable[source] = true;
    reachable[target] = true;

    return reachable;
}

}  // namespace rcspp::preprocess
