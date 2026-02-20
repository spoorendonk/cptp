#pragma once

#include <limits>
#include <queue>
#include <vector>

#include "core/problem.h"

namespace cptp::preprocess {

/// Dijkstra from depot with weight = demand(target_node).
/// Returns reachable[v] = true iff round-trip demand bound is feasible:
///   2 * min_demand_path(depot, v) - demand(v) <= Q
inline std::vector<bool> demand_reachability(const Problem& prob) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();

    constexpr double inf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, inf);
    dist[depot] = 0.0;

    // (cumulative_demand, node)
    using Entry = std::pair<double, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
    pq.emplace(0.0, depot);

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

    // Round-trip check: 2*dist[v] - demand(v) <= Q
    std::vector<bool> reachable(n, false);
    for (int32_t v = 0; v < n; ++v) {
        if (dist[v] < inf) {
            reachable[v] = (2.0 * dist[v] - prob.demand(v) <= Q);
        }
    }
    // Depot is always reachable
    reachable[depot] = true;

    return reachable;
}

}  // namespace cptp::preprocess
