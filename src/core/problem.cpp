#include "core/problem.h"

#include <stdexcept>
#include <utility>

namespace cptp {

void Problem::build(int32_t num_nodes, std::span<const Edge> edges,
                    std::span<const double> edge_costs, std::span<const double> profits,
                    std::span<const double> demands, double capacity, int32_t source,
                    int32_t target) {
    std::vector<std::vector<double>> rows;
    rows.emplace_back(demands.begin(), demands.end());
    build(num_nodes, edges, std::vector<double>(edge_costs.begin(), edge_costs.end()),
          std::vector<double>(profits.begin(), profits.end()), std::move(rows), {capacity}, source,
          target);
}

void Problem::build(int32_t num_nodes, std::span<const Edge> edges, std::vector<double> edge_costs,
                    std::vector<double> profits, std::vector<double> demands, double capacity,
                    int32_t source, int32_t target) {
    std::vector<std::vector<double>> rows;
    rows.push_back(std::move(demands));
    build(num_nodes, edges, std::move(edge_costs), std::move(profits), std::move(rows), {capacity},
          source, target);
}

void Problem::build(int32_t num_nodes, std::span<const Edge> edges, std::vector<double> edge_costs,
                    std::vector<double> profits, std::vector<std::vector<double>> demand_rows,
                    std::vector<double> capacities, int32_t source, int32_t target) {
    if (profits.size() != static_cast<size_t>(num_nodes)) {
        throw std::invalid_argument("profits size must equal num_nodes");
    }
    if (demand_rows.empty()) {
        throw std::invalid_argument("at least one demand row is required");
    }
    if (capacities.size() != demand_rows.size()) {
        throw std::invalid_argument("one capacity per demand row is required");
    }
    for (const auto& row : demand_rows) {
        if (row.size() != static_cast<size_t>(num_nodes)) {
            throw std::invalid_argument("demands size must equal num_nodes");
        }
    }
    if (edge_costs.size() != edges.size()) {
        throw std::invalid_argument("edge_costs size must equal edges size");
    }

    std::vector<static_graph::Edge> graph_edges(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        int32_t u = edges[i].tail;
        int32_t v = edges[i].head;
        if (u > v) {
            std::swap(u, v);
        }
        if (u == v) {
            throw std::invalid_argument("self-loop edge not allowed");
        }
        graph_edges[i] = {u, v};
    }

    num_nodes_ = num_nodes;
    source_ = source;
    target_ = target;
    capacities_ = std::move(capacities);
    graph_ = static_graph(num_nodes, graph_edges);
    edge_costs_ = std::move(edge_costs);
    num_edges_ = static_cast<int32_t>(edges.size());
    profits_ = std::move(profits);
    demands_ = std::move(demand_rows);
}

bool Problem::fits_alone(int32_t v) const {
    for (size_t k = 0; k < capacities_.size(); ++k) {
        if (demands_[k][v] > capacities_[k]) {
            return false;
        }
    }
    return true;
}

}  // namespace cptp
