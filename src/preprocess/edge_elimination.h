#pragma once

#include "core/problem.h"
#include "sep/capacity_rows.h"
#include "sep/separator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cptp::preprocess {

/// No predecessor: the root label, which may extend to any neighbour.
inline constexpr int32_t kNoPred = -1;

/// Label for the capacity-aware labeling.
/// Tracks net cost (edge costs - collected profits) and accumulated demand.
struct Label {
    double cost;    // net cost = sum(edge_costs) - sum(profits along path)
    double demand;  // accumulated demand along path
};

/// Each node's final Pareto front of a labeling, as (cost, demand) labels:
/// every walk the labeling relaxes reaches the node at a cost and a demand
/// no lower than some label of it (#144).
using LabelFronts = std::vector<std::vector<Label>>;

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

/// Parse the `labeling_join` option (#144): "min" (each labeling half at its
/// cheapest label) or "load" (the halves joined by load, edge_lower_bounds);
/// true for "load".
inline bool parse_join_by_load(std::string_view value) {
    if (value == "min" || value == "load") {
        return value == "load";
    }
    throw std::invalid_argument("labeling_join must be min or load, got '" + std::string(value) +
                                "'");
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

/// When a labeling stops: the solve's time limit (#115, #143). The pop budget
/// alone does not bound a labeling's work, since each pop scans every live
/// label at each neighbour, so it can run far past time_limit.
using LabelingDeadline = std::chrono::steady_clock::time_point;
inline constexpr LabelingDeadline kNoLabelingDeadline = LabelingDeadline::max();

/// Whether a labeling has passed its deadline. The clock is read once per
/// kCheckWork labels scanned, the labeling's real work, so the check costs
/// little and still fires within a fraction of a second when one pop scans
/// long live lists. It sets *timed_out (when given) once it fires.
class DeadlineCheck {
public:
    DeadlineCheck(LabelingDeadline deadline, bool* timed_out)
        : deadline_(deadline), timed_out_(timed_out) {}

    /// Counts `scanned` labels of work; true once the deadline has passed.
    bool passed(int64_t scanned) {
        if (deadline_ == kNoLabelingDeadline) {
            return false;
        }
        work_ += scanned;
        if (work_ < next_check_) {
            return false;
        }
        next_check_ = work_ + kCheckWork;
        if (std::chrono::steady_clock::now() < deadline_) {
            return false;
        }
        if (timed_out_ != nullptr) {
            *timed_out_ = true;
        }
        return true;
    }

private:
    static constexpr int64_t kCheckWork = int64_t{1} << 14;
    LabelingDeadline deadline_;
    bool* timed_out_;
    int64_t work_ = 0;
    int64_t next_check_ = 0;
};

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
/// budget of max_queue_pops (0 = none) runs out, when the deadline passes
/// (*timed_out set), or up front when the row has a zero-demand negative
/// cycle (has_zero_demand_negative_cycle), on which the labeling would not
/// terminate. An early stop gives no bound rather than a partial one, whose
/// f[v] would be too high.
///
/// With two_cycle (#144) no walk goes straight back to the node it came from
/// (u -> v -> u), so f[v] bounds the 2-cycle-free walks, which still include
/// every elementary path. A label keeps its predecessor; one label dominates
/// another at the same node when it is <= in cost and demand and has the same
/// predecessor (or none, the root's), and two <= labels with different
/// predecessors dominate together, since one of them can always go where the
/// dominated one would. Only the one-label rule erases a live label.
///
/// When `fronts` is given, (*fronts)[v] is node v's final Pareto front (its
/// live labels, LabelFronts) after a complete run, and empty after an early
/// stop.
inline std::vector<double> labeling_demand(
    const Problem& prob, int32_t root, std::span<const double> demand, double Q,
    std::span<const double> edge_costs, std::span<const double> profits, int64_t max_queue_pops = 0,
    LabelingDeadline deadline = kNoLabelingDeadline, bool* timed_out = nullptr,
    bool two_cycle = false, LabelFronts* fronts = nullptr) {
    if (fronts != nullptr) {
        fronts->clear();
    }
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
    std::vector<int32_t> pred_of;  // each label's predecessor, kept only under two_cycle
    std::vector<char> dead;
    std::vector<std::vector<int32_t>> live(static_cast<size_t>(n));
    std::vector<double> best_cost(static_cast<size_t>(n), kInf);
    std::vector<int32_t> queue;
    auto pred = [&](int32_t id) { return pred_of[static_cast<size_t>(id)]; };

    auto push = [&](int32_t v, Label lbl, int32_t from) {
        const auto id = static_cast<int32_t>(pool.size());
        pool.push_back(lbl);
        node_of.push_back(v);
        if (two_cycle) {
            pred_of.push_back(from);
        }
        dead.push_back(0);
        live[static_cast<size_t>(v)].push_back(id);
        best_cost[static_cast<size_t>(v)] = std::min(lbl.cost, best_cost[static_cast<size_t>(v)]);
        queue.push_back(id);
    };

    push(root, {-profits[root], demand[root]}, kNoPred);
    DeadlineCheck out_of_time(deadline, timed_out);

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
        const int32_t lbl_pred = two_cycle ? pred(id) : kNoPred;

        for (auto e : graph.incident_edges(u)) {
            const int32_t v = graph.other_endpoint(e, u);
            if (two_cycle && v == lbl_pred) {
                continue;  // u -> v -> u
            }

            const double new_demand = lbl.demand + demand[v];
            if (new_demand > Q) {
                continue;
            }

            const double new_cost = lbl.cost + edge_costs[e] - profits[v];

            auto& at_v = live[static_cast<size_t>(v)];
            if (out_of_time.passed(static_cast<int64_t>(at_v.size()) + 1)) {
                return {};  // out of time — no bound, as for the budget
            }
            bool dominated = false;
            int32_t other_pred = kNoPred;  // pred of a <= label that cannot dominate alone
            for (const int32_t other : at_v) {
                const Label& existing = pool[static_cast<size_t>(other)];
                if (existing.cost > new_cost || existing.demand > new_demand) {
                    continue;
                }
                if (!two_cycle || pred(other) == u || pred(other) == kNoPred) {
                    dominated = true;
                    break;
                }
                if (other_pred == kNoPred) {
                    other_pred = pred(other);
                } else if (pred(other) != other_pred) {
                    dominated = true;  // two preds: one of them goes wherever the new one would
                    break;
                }
            }
            if (dominated) {
                continue;
            }

            std::erase_if(at_v, [&](int32_t other) {
                const Label& ex = pool[static_cast<size_t>(other)];
                if (new_cost <= ex.cost && new_demand <= ex.demand &&
                    (!two_cycle || pred(other) == u)) {
                    dead[static_cast<size_t>(other)] = 1;
                    return true;
                }
                return false;
            });

            push(v, {new_cost, new_demand}, u);
        }
    }

    if (fronts != nullptr) {
        fronts->resize(static_cast<size_t>(n));
        for (size_t v = 0; v < live.size(); ++v) {
            auto& front = (*fronts)[v];
            front.reserve(live[v].size());
            for (const int32_t id : live[v]) {
                front.push_back(pool[static_cast<size_t>(id)]);
            }
        }
    }
    return best_cost;
}

/// labeling_demand under capacity row k alone.
inline std::vector<double> labeling_row(const Problem& prob, int32_t root, int32_t k,
                                        std::span<const double> edge_costs,
                                        std::span<const double> profits, int64_t max_queue_pops = 0,
                                        LabelingDeadline deadline = kNoLabelingDeadline,
                                        bool* timed_out = nullptr, bool two_cycle = false,
                                        LabelFronts* fronts = nullptr) {
    return labeling_demand(prob, root, prob.demands(k), prob.capacity(k), edge_costs, profits,
                           max_queue_pops, deadline, timed_out, two_cycle, fronts);
}

/// The surrogate capacity row (#102), shared with the capacity separators
/// (#112): see sep::surrogate_row.
using sep::surrogate_row;

/// Labels the vector labeling (labeling_vector) created, and whether any of
/// labeling_from's labelings stopped at its deadline, for the log.
struct LabelingStats {
    int64_t labels = 0;
    bool timed_out = false;
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
/// Empty when the budget of max_queue_pops (0 = none) runs out, when the
/// deadline passes (stats->timed_out set), when no row is bounded, or when
/// the nodes with zero demand on every bounded row carry a negative cycle
/// (has_zero_demand_negative_cycle on the row sums): only those cycles add
/// no load. An early stop gives no bound rather than a partial one, whose
/// f[v] would be too high.
inline std::vector<double> labeling_vector(const Problem& prob, int32_t root,
                                           std::span<const double> edge_costs,
                                           std::span<const double> profits,
                                           int64_t max_queue_pops = 0,
                                           LabelingStats* stats = nullptr,
                                           LabelingDeadline deadline = kNoLabelingDeadline) {
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
    DeadlineCheck out_of_time(deadline, stats != nullptr ? &stats->timed_out : nullptr);

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
            if (out_of_time.passed(static_cast<int64_t>(at_v.size()) + 1)) {
                return {};  // out of time — no bound, as for the budget
            }
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
/// label count goes to `stats` when given. Every labeling stops at
/// `deadline` (#143): row 0 then gives no bound, so neither does
/// labeling_from, and a later row is skipped (row 0 alone then);
/// stats->timed_out records either. two_cycle (#144) applies to every
/// single-row labeling (labeling_demand), not to the vector one, and
/// `fronts` gets row 0's fronts.
inline std::vector<double> labeling_from(
    const Problem& prob, int32_t root, std::span<const double> edge_costs,
    std::span<const double> profits, int64_t max_queue_pops = 0, LabelRows rows = LabelRows::kFirst,
    int64_t later_rows_max_queue_pops = -1, LabelingStats* stats = nullptr,
    LabelingDeadline deadline = kNoLabelingDeadline, bool two_cycle = false,
    LabelFronts* fronts = nullptr) {
    bool* const timed_out = stats != nullptr ? &stats->timed_out : nullptr;
    auto best = labeling_row(prob, root, 0, edge_costs, profits, max_queue_pops, deadline,
                             timed_out, two_cycle, fronts);
    if (rows == LabelRows::kFirst || best.empty()) {
        return best;
    }
    const int64_t later_pops =
        later_rows_max_queue_pops < 0 ? max_queue_pops : later_rows_max_queue_pops;
    if (rows == LabelRows::kVector) {
        const auto row =
            labeling_vector(prob, root, edge_costs, profits, later_pops, stats, deadline);
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
        const auto row =
            labeling_demand(prob, root, surrogate.demand, surrogate.capacity, edge_costs, profits,
                            later_pops, deadline, timed_out, two_cycle);
        for (size_t v = 0; v < row.size(); ++v) {  // empty row: row 0 alone
            best[v] = std::max(best[v], row[v]);
        }
        return best;
    }
    for (int32_t k = 1; k < prob.num_capacities(); ++k) {
        if (prob.capacity(k) >= sep::kInfiniteCapacity) {
            continue;  // no limit: its labeling is never above row 0's
        }
        const auto row = labeling_row(prob, root, k, edge_costs, profits, later_pops, deadline,
                                      timed_out, two_cycle);
        if (row.empty()) {
            continue;  // budget, deadline or zero-demand negative cycle: no bound from row k
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
                                         LabelingStats* stats = nullptr,
                                         LabelingDeadline deadline = kNoLabelingDeadline,
                                         bool two_cycle = false, LabelFronts* fronts = nullptr) {
    return labeling_from(prob, root, prob.edge_costs(), prob.profits(), max_queue_pops, rows,
                         later_rows_max_queue_pops, stats, deadline, two_cycle, fronts);
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

/// A node's front as a staircase: demand ascending, cost strictly
/// descending. Each kept label is the cheapest of the front at its demand or
/// below, so the cheapest label within a demand budget is the last that fits.
inline std::vector<Label> staircase(std::vector<Label> front) {
    std::ranges::sort(front, [](const Label& a, const Label& b) {
        return a.demand < b.demand || (a.demand == b.demand && a.cost < b.cost);
    });
    std::vector<Label> out;
    for (const Label& l : front) {
        if (out.empty() || l.cost < out.back().cost) {
            out.push_back(l);
        }
    }
    return out;
}

/// Cheapest a.cost + b.cost over labels a of `a`, b of `b` (staircases) with
/// a.demand + b.demand <= cap; +inf when no pair fits.
inline double join_staircases(std::span<const Label> a, std::span<const Label> b, double cap) {
    double best = std::numeric_limits<double>::infinity();
    size_t j = b.size();
    for (const Label& x : a) {  // demand ascending, so the budget left for b shrinks
        while (j > 0 && x.demand + b[j - 1].demand > cap) {
            --j;
        }
        if (j == 0) {
            break;
        }
        best = std::min(best, x.cost + b[j - 1].cost);
    }
    return best;
}

/// Lower bound lb_e on the net cost of any tour or path using edge e, from
/// the labeling bounds f (source -> v) and b (v -> target): the better
/// orientation of f[u] + c_e + b[v] + correction, +inf when neither side
/// reaches (correction as edge_elimination).
///
/// With row 0's fronts from source and target (front_f, front_b; non-null
/// and non-empty, LabelFronts), the two halves are also joined by load
/// (#144): a path source -> u and one v -> target that together form a
/// tour or path carry at most Q_0 of row 0, plus the depot's demand, which
/// both halves of a tour count. So
///   lb_e >= min over a in F_s(u), b in F_t(v), a.demand + b.demand <= cap
///           of a.cost + c_e + b.cost + correction,
/// with cap = Q_0 + d_0(depot) for a tour and Q_0 for a path, again over
/// both orientations. lb_e is the larger of the two bounds: f and b may
/// come from more rows than row 0 (labeling_from). The join stops at
/// `deadline` (the solve's time limit, as the labelings): the edges it has
/// not reached keep the first bound, which is still valid.
inline std::vector<double> edge_lower_bounds(const Problem& prob, std::span<const double> f,
                                             std::span<const double> b, double correction,
                                             const LabelFronts* front_f = nullptr,
                                             const LabelFronts* front_b = nullptr,
                                             LabelingDeadline deadline = kNoLabelingDeadline) {
    const auto& graph = prob.graph();
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::vector<double> lb(static_cast<size_t>(prob.num_edges()), kInf);
    for (auto e : graph.edges()) {
        const int32_t u = graph.edge_source(e);
        const int32_t v = graph.edge_target(e);
        const double c = prob.edge_cost(e);
        lb[static_cast<size_t>(e)] = std::min(f[u] + c + b[v], f[v] + c + b[u]) + correction;
    }
    if (front_f == nullptr || front_b == nullptr || front_f->empty() || front_b->empty()) {
        return lb;
    }
    // A little slack: the halves' loads are sums in a different order than
    // the tour's, so a tour at exactly Q_0 must not be cut off by rounding.
    const double exact_cap =
        prob.capacity(0) + (prob.is_tour() ? prob.demand(0, prob.source()) : 0.0);
    const double cap = exact_cap + 1e-9 * std::max(1.0, std::abs(exact_cap));
    auto stairs = [](const LabelFronts& fronts) {
        std::vector<std::vector<Label>> out;
        out.reserve(fronts.size());
        for (const auto& front : fronts) {
            out.push_back(staircase(front));
        }
        return out;
    };
    const auto sf = stairs(*front_f);
    const auto own_sb = front_f == front_b ? std::vector<std::vector<Label>>{} : stairs(*front_b);
    const auto& sb = front_f == front_b ? sf : own_sb;  // a tour shares one set
    DeadlineCheck out_of_time(deadline, nullptr);
    for (auto e : graph.edges()) {
        const int32_t u = graph.edge_source(e);
        const int32_t v = graph.edge_target(e);
        const size_t work = sf[u].size() + sb[v].size() + sf[v].size() + sb[u].size();
        if (out_of_time.passed(static_cast<int64_t>(work))) {
            break;  // out of time: the rest keep the first bound
        }
        const double c = prob.edge_cost(e);
        const double joined =
            std::min(join_staircases(sf[u], sb[v], cap), join_staircases(sf[v], sb[u], cap)) + c +
            correction;
        auto& x = lb[static_cast<size_t>(e)];
        x = std::max(x, joined);
    }
    return lb;
}

/// Raises edge_lb (edge_lower_bounds) to the plain bound from f and b where
/// that is higher, after f and b were recomputed (all-pairs bounds).
inline void raise_edge_lower_bounds(std::vector<double>& edge_lb, const Problem& prob,
                                    std::span<const double> f, std::span<const double> b,
                                    double correction) {
    const auto plain = edge_lower_bounds(prob, f, b, correction);
    for (size_t e = 0; e < edge_lb.size(); ++e) {
        edge_lb[e] = std::max(edge_lb[e], plain[e]);
    }
}

/// Edge elimination from per-edge lower bounds (edge_lower_bounds):
/// eliminated[e] when lb_e > upper_bound, none when the bound is +inf.
inline std::vector<bool> edge_elimination(const Problem& prob, std::span<const double> edge_lb,
                                          double upper_bound) {
    std::vector<bool> eliminated(static_cast<size_t>(prob.num_edges()), false);
    if (upper_bound >= std::numeric_limits<double>::infinity()) {
        return eliminated;
    }
    for (size_t e = 0; e < eliminated.size(); ++e) {
        eliminated[e] = edge_lb[e] > upper_bound + 1e-6;
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
