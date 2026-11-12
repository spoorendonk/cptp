#pragma once

#include "core/problem.h"

#include <limits>
#include <queue>
#include <vector>

namespace cptp::preprocess {

/// Dijkstra from a given root with weight = demand(k, target_node) of
/// capacity row k. Returns dist[v] = minimum cumulative row-k demand on any
/// path from root to v.
inline std::vector<double> demand_dijkstra(const Problem& prob, int32_t root, int32_t k = 0) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    std::vector<double> dist(n, kInf);
    dist[root] = 0.0;

    using Entry = std::pair<double, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
    pq.emplace(0.0, root);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) {
            continue;
        }

        for (auto e : graph.incident_edges(u)) {
            int32_t v = graph.other_endpoint(e, u);
            double nd = d + prob.demand(k, v);
            if (nd < dist[v]) {
                dist[v] = nd;
                pq.emplace(nd, v);
            }
        }
    }

    return dist;
}

/// Demand-reachability preprocessing for capacity row k.
/// Tour (source == target): round-trip check 2*dist[v] - demand(v) <= Q.
/// Path (source != target): one-way check dist_s[v] + dist_t[v] - demand(v) <=
/// Q. A node unreachable under any one row is unreachable under all of them.
inline std::vector<bool> demand_reachability(const Problem& prob, int32_t k = 0) {
    const int32_t n = prob.num_nodes();
    const double Q = prob.capacity(k);
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    auto dist_s = demand_dijkstra(prob, source, k);

    std::vector<bool> reachable(n, false);

    if (prob.is_tour()) {
        // Tour: round-trip check
        for (int32_t v = 0; v < n; ++v) {
            if (dist_s[v] < kInf) {
                reachable[v] = (2.0 * dist_s[v] - prob.demand(k, v) <= Q);
            }
        }
    } else {
        // Path: Dijkstra from both source and target
        auto dist_t = demand_dijkstra(prob, target, k);
        for (int32_t v = 0; v < n; ++v) {
            if (dist_s[v] < kInf && dist_t[v] < kInf) {
                reachable[v] = (dist_s[v] + dist_t[v] - prob.demand(k, v) <= Q);
            }
        }
    }

    // Source and target are always reachable
    reachable[source] = true;
    reachable[target] = true;

    return reachable;
}

}  // namespace cptp::preprocess
