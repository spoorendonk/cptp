#pragma once

#include "core/problem.h"
#include "sep/capacity_rows.h"
#include "sep/separator.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cptp::preprocess {

/// Label for the capacity-aware labeling.
/// Tracks net cost (edge costs - collected profits) and accumulated demand.
struct Label {
    double cost;    // net cost = sum(edge_costs) - sum(profits along path)
    double demand;  // accumulated demand along path
};

/// Which capacity rows the labeling bounds come from (#101).
enum class LabelRows {
    kFirst,      ///< row 0 only (the single-row behaviour)
    kAll,        ///< row 0 and every row k >= 1 with Q_k < kInfiniteCapacity, combined
    kSurrogate,  ///< row 0 combined with one surrogate row (#102, surrogate_row)
    kVector,     ///< row 0 combined with vector labels over every bounded row (#103)
};

/// Parse the `labeling_rows` option: "first", "all", "surrogate" or "vector".
inline LabelRows parse_label_rows(std::string_view value) {
    if (value == "first") {
        return LabelRows::kFirst;
    }
    if (value == "all") {
        return LabelRows::kAll;
    }
    if (value == "surrogate") {
        return LabelRows::kSurrogate;
    }
    if (value == "vector") {
        return LabelRows::kVector;
    }
    throw std::invalid_argument("labeling_rows must be first, all, surrogate or vector, got '" +
                                std::string(value) + "'");
}

/// Whether the zero-demand nodes of the row `demand` carry a closed walk of negative net cost
/// (arc u -> v costs edge_costs[e] - profits[v]). Such a walk adds no demand,
/// so the non-elementary labeling under that row would improve its labels along
/// it forever. Without one the labeling terminates: every other cycle adds
/// demand, and a zero-cost one is caught by the dominance test.
///
/// The graph is undirected, so a closed walk enters each node as often as it
/// leaves it, and its cost is sum over its arcs of c_e - (p_u + p_v) / 2: half
/// the sum of its edges' 2-cycle costs 2 c_e - p_u - p_v. A negative closed
/// walk exists exactly when some zero-demand edge's 2-cycle is negative.
inline bool has_zero_demand_negative_cycle(const Problem& prob, std::span<const double> demand,
                                           std::span<const double> edge_costs,
                                           std::span<const double> profits) {
    const auto& graph = prob.graph();
    for (int32_t e = 0; e < graph.num_edges(); ++e) {
        const int32_t u = graph.edge_source(e);
        const int32_t v = graph.edge_target(e);
        if (demand[u] <= 0.0 && demand[v] <= 0.0 &&
            2.0 * edge_costs[e] - profits[u] - profits[v] < 0.0) {
            return true;
        }
    }
    return false;
}

/// has_zero_demand_negative_cycle for capacity row k.
inline bool has_zero_demand_negative_cycle(const Problem& prob, int32_t k,
                                           std::span<const double> edge_costs,
                                           std::span<const double> profits) {
    return has_zero_demand_negative_cycle(prob, prob.demands(k), edge_costs, profits);
}

/// Capacity-aware labeling from a given root node under one capacity row,
/// demand[v] against Q, with custom edge costs and node profits.
///
/// For each node v, computes the minimum net cost f[v] to reach v from root
/// via any (possibly non-elementary) path that respects that row's capacity.
///
/// Net cost = sum of edge costs along path - sum of profits of visited nodes
/// (excluding root profit, which is always collected... see seed below).
///
/// Returns f[v] for each node, f[root] = -profits[root]; empty when the
/// budget of max_queue_pops (0 = none) runs out, or up front when the row has
/// a zero-demand negative cycle (has_zero_demand_negative_cycle), on which
/// the labeling would not terminate.
inline std::vector<double> labeling_demand(const Problem& prob, int32_t root,
                                           std::span<const double> demand, double Q,
                                           std::span<const double> edge_costs,
                                           std::span<const double> profits,
                                           int64_t max_queue_pops = 0) {
    if (has_zero_demand_negative_cycle(prob, demand, edge_costs, profits)) {
        return {};
    }
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    // Labels live in one pool under stable ids, a dominated one marked dead,
    // and the queue holds ids, as in labeling_vector: queueing (node, index)
    // would let erase_if shift a queued label's index onto another, leaving a
    // live label unextended and f too high (#111).
    std::vector<Label> pool;
    std::vector<int32_t> node_of;
    std::vector<char> dead;
    std::vector<std::vector<int32_t>> live(static_cast<size_t>(n));
    std::vector<double> best_cost(static_cast<size_t>(n), kInf);
    std::vector<int32_t> queue;

    auto push = [&](int32_t v, Label lbl) {
        const auto id = static_cast<int32_t>(pool.size());
        pool.push_back(lbl);
        node_of.push_back(v);
        dead.push_back(0);
        live[static_cast<size_t>(v)].push_back(id);
        best_cost[static_cast<size_t>(v)] = std::min(lbl.cost, best_cost[static_cast<size_t>(v)]);
        queue.push_back(id);
    };

    push(root, {-profits[root], demand[root]});

    size_t head = 0;
    while (head < queue.size()) {
        if (max_queue_pops > 0 && static_cast<int64_t>(head) >= max_queue_pops) {
            return {};  // budget exhausted — caller checks .empty()
        }
        const int32_t id = queue[head++];
        if (dead[static_cast<size_t>(id)] != 0) {
            continue;
        }
        const int32_t u = node_of[static_cast<size_t>(id)];
        const Label lbl = pool[static_cast<size_t>(id)];  // a copy: push grows the pool

        for (auto e : graph.incident_edges(u)) {
            const int32_t v = graph.other_endpoint(e, u);

            const double new_demand = lbl.demand + demand[v];
            if (new_demand > Q) {
                continue;
            }

            const double new_cost = lbl.cost + edge_costs[e] - profits[v];

            auto& at_v = live[static_cast<size_t>(v)];
            bool dominated = false;
            for (const int32_t other : at_v) {
                const Label& existing = pool[static_cast<size_t>(other)];
                if (existing.cost <= new_cost && existing.demand <= new_demand) {
                    dominated = true;
                    break;
                }
            }
            if (dominated) {
                continue;
            }

            std::erase_if(at_v, [&](int32_t other) {
                const Label& ex = pool[static_cast<size_t>(other)];
                if (new_cost <= ex.cost && new_demand <= ex.demand) {
                    dead[static_cast<size_t>(other)] = 1;
                    return true;
                }
                return false;
            });

            push(v, {new_cost, new_demand});
        }
    }

    return best_cost;
}

/// labeling_demand under capacity row k alone.
inline std::vector<double> labeling_row(const Problem& prob, int32_t root, int32_t k,
                                        std::span<const double> edge_costs,
                                        std::span<const double> profits,
                                        int64_t max_queue_pops = 0) {
    return labeling_demand(prob, root, prob.demands(k), prob.capacity(k), edge_costs, profits,
                           max_queue_pops);
}

/// The surrogate capacity row (#102), shared with the capacity separators
/// (#112): see sep::surrogate_row.
using sep::surrogate_row;

/// Labels the vector labeling (labeling_vector) created, for the log.
struct LabelingStats {
    int64_t labels = 0;
};

/// Capacity-aware labeling from root with vector labels: each label carries
/// its load on every row with Q_k < kInfiniteCapacity, a label extends only
/// when every load stays within its Q_k, and one label dominates another
/// when its cost and every load are <=. So f[v] is the minimum net cost of
/// a (possibly non-elementary) path that respects all rows at once, a bound
/// at least as strong as each row's own and the surrogate's.
///
/// Labels live in one pool under stable ids, a dominated one marked dead,
/// and the queue holds ids, so every label live when popped is extended.
/// Empty when the budget of max_queue_pops (0 = none) runs out, when no row
/// is bounded, or when the nodes with zero demand on every bounded row carry
/// a negative cycle (has_zero_demand_negative_cycle on the row sums): only
/// those cycles add no load.
inline std::vector<double> labeling_vector(const Problem& prob, int32_t root,
                                           std::span<const double> edge_costs,
                                           std::span<const double> profits,
                                           int64_t max_queue_pops = 0,
                                           LabelingStats* stats = nullptr) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    std::vector<int32_t> rows;
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        if (prob.capacity(k) < sep::kInfiniteCapacity) {
            rows.push_back(k);
        }
    }
    if (rows.empty()) {
        return {};
    }
    std::vector<double> total(static_cast<size_t>(n), 0.0);
    for (const int32_t k : rows) {
        for (int32_t v = 0; v < n; ++v) {
            total[static_cast<size_t>(v)] += prob.demand(k, v);
        }
    }
    if (has_zero_demand_negative_cycle(prob, total, edge_costs, profits)) {
        return {};
    }

    const size_t r = rows.size();
    const size_t stride = r + 1;  // pool[id * stride] = cost, then r loads
    std::vector<double> pool;
    std::vector<int32_t> node_of;
    std::vector<char> dead;
    std::vector<std::vector<int32_t>> live(static_cast<size_t>(n));
    std::vector<double> best_cost(static_cast<size_t>(n), kInf);
    std::vector<int32_t> queue;

    auto push = [&](int32_t v, double cost, std::span<const double> loads) {
        const auto id = static_cast<int32_t>(node_of.size());
        pool.push_back(cost);
        pool.insert(pool.end(), loads.begin(), loads.end());
        node_of.push_back(v);
        dead.push_back(0);
        live[static_cast<size_t>(v)].push_back(id);
        best_cost[static_cast<size_t>(v)] = std::min(cost, best_cost[static_cast<size_t>(v)]);
        queue.push_back(id);
        if (stats != nullptr) {
            ++stats->labels;
        }
    };

    std::vector<double> loads(r);
    for (size_t j = 0; j < r; ++j) {
        loads[j] = prob.demand(rows[j], root);
    }
    push(root, -profits[root], loads);

    size_t head = 0;
    while (head < queue.size()) {
        if (max_queue_pops > 0 && static_cast<int64_t>(head) >= max_queue_pops) {
            return {};  // budget exhausted — caller checks .empty()
        }
        const int32_t id = queue[head++];
        if (dead[static_cast<size_t>(id)] != 0) {
            continue;
        }
        const int32_t u = node_of[static_cast<size_t>(id)];
        const size_t base = static_cast<size_t>(id) * stride;
        const double cost = pool[base];
        std::vector<double> from(pool.begin() + static_cast<ptrdiff_t>(base) + 1,
                                 pool.begin() + static_cast<ptrdiff_t>(base + stride));

        for (auto e : graph.incident_edges(u)) {
            const int32_t v = graph.other_endpoint(e, u);
            bool fits = true;
            for (size_t j = 0; j < r; ++j) {
                loads[j] = from[j] + prob.demand(rows[j], v);
                if (loads[j] > prob.capacity(rows[j])) {
                    fits = false;
                    break;
                }
            }
            if (!fits) {
                continue;
            }
            const double new_cost = cost + edge_costs[e] - profits[v];

            // a <= b componentwise, over cost and every load
            auto leq = [&](double cost_a, const double* a, double cost_b, const double* b) {
                if (cost_a > cost_b) {
                    return false;
                }
                for (size_t j = 0; j < r; ++j) {
                    if (a[j] > b[j]) {
                        return false;
                    }
                }
                return true;
            };
            auto& at_v = live[static_cast<size_t>(v)];
            bool dominated = false;
            for (const int32_t other : at_v) {
                const size_t ob = static_cast<size_t>(other) * stride;
                if (leq(pool[ob], &pool[ob + 1], new_cost, loads.data())) {
                    dominated = true;
                    break;
                }
            }
            if (dominated) {
                continue;
            }
            std::erase_if(at_v, [&](int32_t other) {
                const size_t ob = static_cast<size_t>(other) * stride;
                if (leq(new_cost, loads.data(), pool[ob], &pool[ob + 1])) {
                    dead[static_cast<size_t>(other)] = 1;
                    return true;
                }
                return false;
            });
            push(v, new_cost, loads);
        }
    }

    return best_cost;
}

/// Capacity-aware labeling from root (see labeling_row) over the rows `rows`
/// names. Each row's f_k[v] is a lower bound on the net cost of reaching v
/// under all rows (it relaxes the others), so their maximum is one too:
/// f[v] = max_k f_k[v]. Row 0 alone is the single-row labeling. Empty when
/// row 0's labeling comes back empty (labeling_row); a row k >= 1 that does
/// is skipped, since the maximum over the other rows is still a bound, so
/// kAll is never weaker than kFirst. Rows k >= 1 get their own budget,
/// later_rows_max_queue_pops (-1 = max_queue_pops), so a caller that runs row
/// 0 unbounded need not run the others unbounded too. kSurrogate combines
/// row 0 with the surrogate row (surrogate_row) the same way, under that
/// budget, and kVector with the vector labeling (labeling_vector), whose
/// label count goes to `stats` when given.
inline std::vector<double> labeling_from(
    const Problem& prob, int32_t root, std::span<const double> edge_costs,
    std::span<const double> profits, int64_t max_queue_pops = 0, LabelRows rows = LabelRows::kFirst,
    int64_t later_rows_max_queue_pops = -1, LabelingStats* stats = nullptr) {
    auto best = labeling_row(prob, root, 0, edge_costs, profits, max_queue_pops);
    if (rows == LabelRows::kFirst || best.empty()) {
        return best;
    }
    const int64_t later_pops =
        later_rows_max_queue_pops < 0 ? max_queue_pops : later_rows_max_queue_pops;
    if (rows == LabelRows::kVector) {
        const auto row = labeling_vector(prob, root, edge_costs, profits, later_pops, stats);
        for (size_t v = 0; v < row.size(); ++v) {  // empty row: row 0 alone
            best[v] = std::max(best[v], row[v]);
        }
        return best;
    }
    if (rows == LabelRows::kSurrogate) {
        const auto surrogate = surrogate_row(prob);
        if (surrogate.demand.empty()) {
            return best;
        }
        const auto row = labeling_demand(prob, root, surrogate.demand, surrogate.capacity,
                                         edge_costs, profits, later_pops);
        for (size_t v = 0; v < row.size(); ++v) {  // empty row: row 0 alone
            best[v] = std::max(best[v], row[v]);
        }
        return best;
    }
    for (int32_t k = 1; k < prob.num_capacities(); ++k) {
        if (prob.capacity(k) >= sep::kInfiniteCapacity) {
            continue;  // no limit: its labeling is never above row 0's
        }
        const auto row = labeling_row(prob, root, k, edge_costs, profits, later_pops);
        if (row.empty()) {
            continue;  // budget or zero-demand negative cycle: no bound from row k
        }
        for (size_t v = 0; v < best.size(); ++v) {
            best[v] = std::max(best[v], row[v]);
        }
    }
    return best;
}

/// Capacity-aware labeling using the problem's own costs.
inline std::vector<double> labeling_from(const Problem& prob, int32_t root,
                                         int64_t max_queue_pops = 0,
                                         LabelRows rows = LabelRows::kFirst,
                                         int64_t later_rows_max_queue_pops = -1,
                                         LabelingStats* stats = nullptr) {
    return labeling_from(prob, root, prob.edge_costs(), prob.profits(), max_queue_pops, rows,
                         later_rows_max_queue_pops, stats);
}

/// Forward labeling from a given source node (wrapper for labeling_from).
inline std::vector<double> forward_labeling(const Problem& prob, int32_t source,
                                            int64_t max_queue_pops = 0) {
    return labeling_from(prob, source, max_queue_pops);
}

/// Forward labeling from the problem's source node.
inline std::vector<double> forward_labeling(const Problem& prob, int64_t max_queue_pops = 0) {
    return labeling_from(prob, prob.source(), max_queue_pops);
}

/// Backward labeling from a given target node.
/// For undirected graphs, this is identical to forward_labeling(prob, target).
inline std::vector<double> backward_labeling(const Problem& prob, int32_t target,
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
inline std::vector<bool> edge_elimination(const Problem& prob, const std::vector<double>& f,
                                          const std::vector<double>& b, double upper_bound,
                                          double correction) {
    const auto& graph = prob.graph();
    const int32_t m = prob.num_edges();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    std::vector<bool> eliminated(m, false);

    if (upper_bound >= kInf) {
        return eliminated;
    }

    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);

        if (f[u] >= kInf || b[v] >= kInf) {
            // Try both orientations for undirected
            if (f[v] >= kInf || b[u] >= kInf) {
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
inline std::vector<bool> edge_elimination(const Problem& prob, double upper_bound) {
    const auto& graph = prob.graph();
    const int32_t m = prob.num_edges();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    std::vector<bool> eliminated(m, false);

    if (upper_bound >= kInf) {
        return eliminated;
    }

    // Forward labeling from source
    auto f_s = labeling_from(prob, prob.source());
    if (f_s.empty()) {
        return eliminated;  // zero-demand negative cycle: no bounds
    }

    if (prob.is_tour()) {
        // Tour: f_s[u] + c(u,v) + f_s[v] + profit(depot)
        // (undirected: backward = forward, correct for double depot subtraction)
        double depot_profit = prob.profit(prob.source());

        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);

            if (f_s[u] >= kInf || f_s[v] >= kInf) {
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
        if (f_t.empty()) {
            return eliminated;
        }

        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);

            // Try both orientations: source->...->u->v->...->target and reverse
            double lb1 =
                (f_s[u] < kInf && f_t[v] < kInf) ? f_s[u] + prob.edge_cost(e) + f_t[v] : kInf;
            double lb2 =
                (f_s[v] < kInf && f_t[u] < kInf) ? f_s[v] + prob.edge_cost(e) + f_t[u] : kInf;
            double lb = std::min(lb1, lb2);

            if (lb >= kInf || lb > upper_bound + 1e-6) {
                eliminated[e] = true;
            }
        }
    }

    return eliminated;
}

}  // namespace cptp::preprocess
