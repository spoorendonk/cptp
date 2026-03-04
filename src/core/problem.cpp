#include "core/problem.h"

#include <cassert>

namespace cptp {

void Problem::build(int32_t num_nodes,
                    std::span<const Edge> edges,
                    std::span<const double> edge_costs,
                    std::span<const double> profits,
                    std::span<const double> demands,
                    double capacity,
                    int32_t source,
                    int32_t target) {
    assert(profits.size() == static_cast<size_t>(num_nodes));
    assert(demands.size() == static_cast<size_t>(num_nodes));
    assert(edge_costs.size() == edges.size());

    num_nodes_ = num_nodes;
    source_ = source;
    target_ = target;
    capacity_ = capacity;

    // Build static_graph from edge list
    std::vector<static_graph::Edge> graph_edges(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        int32_t u = edges[i].tail;
        int32_t v = edges[i].head;
        assert(u < v);
        graph_edges[i] = {u, v};
    }

    graph_ = static_graph(num_nodes, graph_edges);
    edge_costs_.assign(edge_costs.begin(), edge_costs.end());
    num_edges_ = static_cast<int32_t>(edges.size());

    profits_.assign(profits.begin(), profits.end());
    demands_.assign(demands.begin(), demands.end());
}

void Problem::build(int32_t num_nodes,
                    std::span<const Edge> edges,
                    std::vector<double> edge_costs,
                    std::vector<double> profits,
                    std::vector<double> demands,
                    double capacity,
                    int32_t source,
                    int32_t target) {
    assert(profits.size() == static_cast<size_t>(num_nodes));
    assert(demands.size() == static_cast<size_t>(num_nodes));
    assert(edge_costs.size() == edges.size());

    num_nodes_ = num_nodes;
    source_ = source;
    target_ = target;
    capacity_ = capacity;

    std::vector<static_graph::Edge> graph_edges(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        int32_t u = edges[i].tail;
        int32_t v = edges[i].head;
        assert(u < v);
        graph_edges[i] = {u, v};
    }

    graph_ = static_graph(num_nodes, graph_edges);
    edge_costs_ = std::move(edge_costs);
    num_edges_ = static_cast<int32_t>(edges.size());

    profits_ = std::move(profits);
    demands_ = std::move(demands);
}

}  // namespace cptp
