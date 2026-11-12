#pragma once

#include "core/static_graph.h"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

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

/// CPTP instance: undirected graph with profits, costs, and m >= 1 capacity
/// rows sum_i d^k_i y_i <= Q_k (k = 0..m-1), all of which a solution must
/// satisfy. The single-row accessors (demand(v), demands(), capacity()) are
/// row 0.
class Problem {
public:
    Problem() = default;

    /// Build from raw data (copies from spans), one capacity row.
    /// edges[i] = {tail, head} with tail < head.
    void build(int32_t num_nodes, std::span<const Edge> edges, std::span<const double> edge_costs,
               std::span<const double> profits, std::span<const double> demands, double capacity,
               int32_t source = 0, int32_t target = 0);

    /// Build by moving owned vectors (avoids copies when caller no longer needs
    /// them), one capacity row.
    void build(int32_t num_nodes, std::span<const Edge> edges, std::vector<double> edge_costs,
               std::vector<double> profits, std::vector<double> demands, double capacity,
               int32_t source = 0, int32_t target = 0);

    /// Build with m = demand_rows.size() >= 1 capacity rows: demand_rows[k] has
    /// one entry per node and capacities[k] is its limit.
    void build(int32_t num_nodes, std::span<const Edge> edges, std::vector<double> edge_costs,
               std::vector<double> profits, std::vector<std::vector<double>> demand_rows,
               std::vector<double> capacities, int32_t source = 0, int32_t target = 0);

    // --- Accessors ---
    [[nodiscard]] const Graph& graph() const { return graph_; }
    [[nodiscard]] int32_t num_nodes() const { return num_nodes_; }
    [[nodiscard]] int32_t num_edges() const { return num_edges_; }
    [[nodiscard]] int32_t source() const { return source_; }
    [[nodiscard]] int32_t target() const { return target_; }
    [[nodiscard]] int32_t depot() const { return source_; }  // backward compat alias
    [[nodiscard]] bool is_tour() const { return source_ == target_; }
    /// Number of capacity rows m (>= 1 once built).
    [[nodiscard]] int32_t num_capacities() const {
        return static_cast<int32_t>(capacities_.size());
    }
    [[nodiscard]] double capacity() const { return capacities_[0]; }
    [[nodiscard]] double capacity(int32_t k) const { return capacities_[k]; }
    [[nodiscard]] const std::vector<double>& capacities() const { return capacities_; }

    [[nodiscard]] double edge_cost(int32_t e) const { return edge_costs_[e]; }
    [[nodiscard]] double profit(int32_t v) const { return profits_[v]; }
    [[nodiscard]] double demand(int32_t v) const { return demands_[0][v]; }
    [[nodiscard]] double demand(int32_t k, int32_t v) const { return demands_[k][v]; }

    [[nodiscard]] const std::vector<double>& edge_costs() const { return edge_costs_; }
    [[nodiscard]] const std::vector<double>& profits() const { return profits_; }
    [[nodiscard]] const std::vector<double>& demands() const { return demands_[0]; }
    [[nodiscard]] const std::vector<double>& demands(int32_t k) const { return demands_[k]; }
    [[nodiscard]] const std::vector<std::vector<double>>& demand_rows() const { return demands_; }

    /// True when node v fits alone under every row: d^k_v <= Q_k for all k.
    [[nodiscard]] bool fits_alone(int32_t v) const;

    std::string name;

private:
    Graph graph_;
    int32_t num_nodes_ = 0;
    int32_t num_edges_ = 0;
    int32_t source_ = 0;
    int32_t target_ = 0;
    std::vector<double> capacities_ = {0.0};
    std::vector<double> edge_costs_;
    std::vector<double> profits_;
    std::vector<std::vector<double>> demands_ = {{}};
};

}  // namespace cptp
