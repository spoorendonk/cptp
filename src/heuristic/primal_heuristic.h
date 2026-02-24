#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <span>
#include <thread>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

#include "core/problem.h"

namespace rcspp::heuristic {

/// Subset of the original problem's edges/nodes, built from LP relaxation values.
struct ReducedGraph {
    std::vector<bool> edge_active;  // size m: which edges are in the reduced graph
    std::vector<bool> node_active;  // size n: which nodes are candidates
};

namespace detail {

/// Find edge index between u and v, or -1 if not found.
/// When edge_active is non-empty, inactive edges return -1.
inline int32_t find_edge(const Graph& g, int32_t u, int32_t v,
                         const std::vector<bool>& edge_active = {}) {
    for (int32_t e : g.incident_edges(u)) {
        if (!edge_active.empty() && !edge_active[e]) continue;
        if (g.other_endpoint(e, u) == v) return e;
    }
    return -1;
}

inline double edge_cost(const Problem& prob, int32_t u, int32_t v,
                        const std::vector<bool>& edge_active = {}) {
    int32_t e = find_edge(prob.graph(), u, v, edge_active);
    return (e >= 0) ? prob.edge_cost(e) : std::numeric_limits<double>::max();
}

/// Compute tour/path objective: sum(edge costs) - sum(profits).
/// For tours: subtract depot profit once (depot appears at both ends but counts once).
/// For paths: subtract source and target profits.
inline double tour_objective(const Problem& prob,
                             const std::vector<int32_t>& tour) {
    double cost = 0.0;
    for (size_t i = 0; i + 1 < tour.size(); ++i)
        cost += edge_cost(prob, tour[i], tour[i + 1]);
    for (size_t i = 1; i + 1 < tour.size(); ++i)
        cost -= prob.profit(tour[i]);
    if (prob.is_tour()) {
        // Tour: depot profit counted once (depot is first and last element)
        cost -= prob.profit(prob.source());
    } else {
        // Path: source and target profits (source is first, target is last)
        cost -= prob.profit(tour.front());
        cost -= prob.profit(tour.back());
    }
    return cost;
}

/// Greedy insertion: insert customers in given order at their cheapest position.
/// When edge_active is non-empty, only active edges are considered.
inline void greedy_insert(const Problem& prob,
                          std::vector<int32_t>& tour,
                          std::vector<bool>& in_tour,
                          double& remaining_cap,
                          const std::vector<int32_t>& order,
                          const std::vector<bool>& edge_active = {}) {
    for (int32_t c : order) {
        if (in_tour[c]) continue;
        if (prob.demand(c) > remaining_cap) continue;

        double best_delta = std::numeric_limits<double>::max();
        size_t best_pos = 0;

        for (size_t pos = 1; pos < tour.size(); ++pos) {
            int32_t prev = tour[pos - 1];
            int32_t next = tour[pos];

            double cost_pc = edge_cost(prob, prev, c, edge_active);
            double cost_cn = edge_cost(prob, c, next, edge_active);
            if (cost_pc > 1e15 || cost_cn > 1e15) continue;

            double cost_pn = edge_cost(prob, prev, next, edge_active);
            if (cost_pn > 1e15) cost_pn = 0.0;

            double delta = cost_pc + cost_cn - cost_pn;
            if (delta < best_delta) {
                best_delta = delta;
                best_pos = pos;
            }
        }

        if (best_pos > 0) {
            tour.insert(tour.begin() + static_cast<ptrdiff_t>(best_pos), c);
            in_tour[c] = true;
            remaining_cap -= prob.demand(c);
        }
    }
}

/// Run local search neighborhoods until no improvement or iteration limit.
/// When edge_active is non-empty, only active edges are considered.
inline void local_search(const Problem& prob,
                         std::vector<int32_t>& tour,
                         std::vector<bool>& in_tour,
                         double& remaining_cap,
                         const std::vector<bool>& edge_active = {},
                         int max_iter = 200) {
    const int32_t n = prob.num_nodes();
    // Minimum tour size for node drop: tours need >= 4 (depot + 2 customers + depot),
    // paths need >= 3 (source + 1 customer + target).
    const int32_t min_drop_size = prob.is_tour() ? 5 : 4;

    for (int iter = 0; iter < max_iter; ++iter) {
        bool improved = false;
        int32_t len = static_cast<int32_t>(tour.size());

        // --- 2-opt ---
        for (int32_t i = 0; i < len - 2 && !improved; ++i) {
            for (int32_t j = i + 2; j < len - 1; ++j) {
                double old_c = edge_cost(prob, tour[i], tour[i + 1], edge_active)
                             + edge_cost(prob, tour[j], tour[j + 1], edge_active);
                double new_c = edge_cost(prob, tour[i], tour[j], edge_active)
                             + edge_cost(prob, tour[i + 1], tour[j + 1], edge_active);
                if (new_c < old_c - 1e-9) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improved = true;
                    break;
                }
            }
        }
        if (improved) continue;

        // --- Or-opt (chains of 1, 2, 3) ---
        for (int cl = 1; cl <= 3 && !improved; ++cl) {
            for (int32_t i = 1; i + cl < len - 1 && !improved; ++i) {
                int32_t se = i + cl - 1;
                double remove_save =
                    edge_cost(prob, tour[i - 1], tour[i], edge_active)
                    + edge_cost(prob, tour[se], tour[se + 1], edge_active)
                    - edge_cost(prob, tour[i - 1], tour[se + 1], edge_active);

                for (int32_t j = 1; j < len - 1 && !improved; ++j) {
                    if (j >= i - 1 && j <= se + 1) continue;
                    double insert_add =
                        edge_cost(prob, tour[j - 1], tour[i], edge_active)
                        + edge_cost(prob, tour[se], tour[j], edge_active)
                        - edge_cost(prob, tour[j - 1], tour[j], edge_active);
                    if (insert_add - remove_save < -1e-9) {
                        std::vector<int32_t> seg(tour.begin() + i,
                                                 tour.begin() + se + 1);
                        tour.erase(tour.begin() + i, tour.begin() + se + 1);
                        int32_t ins = (j > se) ? j - cl : j;
                        tour.insert(tour.begin() + ins, seg.begin(), seg.end());
                        improved = true;
                    }
                }
            }
        }
        if (improved) continue;

        // --- Node drop ---
        for (int32_t i = 1; i < len - 1 && len >= min_drop_size && !improved; ++i) {
            int32_t c = tour[i];
            double save = edge_cost(prob, tour[i - 1], c, edge_active)
                        + edge_cost(prob, c, tour[i + 1], edge_active)
                        - edge_cost(prob, tour[i - 1], tour[i + 1], edge_active);
            if (save - prob.profit(c) > 1e-9) {
                tour.erase(tour.begin() + i);
                in_tour[c] = false;
                remaining_cap += prob.demand(c);
                improved = true;
            }
        }
        if (improved) continue;

        // --- Node add ---
        for (int32_t c = 0; c < n && !improved; ++c) {
            if (in_tour[c]) continue;
            if (prob.demand(c) > remaining_cap) continue;

            double best_delta = 1e-9;
            size_t best_pos = 0;

            for (size_t pos = 1; pos < tour.size(); ++pos) {
                double cost_pc = edge_cost(prob, tour[pos - 1], c, edge_active);
                double cost_cn = edge_cost(prob, c, tour[pos], edge_active);
                if (cost_pc > 1e15 || cost_cn > 1e15) continue;
                double cost_pn = edge_cost(prob, tour[pos - 1], tour[pos], edge_active);
                if (cost_pn > 1e15) cost_pn = 0.0;

                double delta = cost_pc + cost_cn - cost_pn - prob.profit(c);
                if (delta < best_delta) {
                    best_delta = delta;
                    best_pos = pos;
                }
            }

            if (best_pos > 0) {
                tour.insert(tour.begin() + static_cast<ptrdiff_t>(best_pos), c);
                in_tour[c] = true;
                remaining_cap -= prob.demand(c);
                improved = true;
            }
        }

        if (!improved) break;
    }
}

/// Single restart: construct + local search, return (tour, objective).
/// When edge_active is non-empty, only active edges are considered.
inline std::pair<std::vector<int32_t>, double>
single_restart(const Problem& prob,
               const std::vector<int32_t>& customers,
               const std::vector<int32_t>& order,
               const std::vector<bool>& edge_active = {}) {
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();
    const bool is_tour = prob.is_tour();

    // Initialize: tour for closed loop, path for open s-t
    std::vector<int32_t> tour;
    std::vector<bool> in_tour(n, false);
    if (is_tour) {
        tour = {source, source};
        in_tour[source] = true;
    } else {
        tour = {source, target};
        in_tour[source] = true;
        in_tour[target] = true;
    }
    // Mandatory endpoints consume capacity too: depot once for tours,
    // source+target for paths.
    double remaining_cap = Q - prob.demand(source);
    if (!is_tour) remaining_cap -= prob.demand(target);

    if (remaining_cap < 0.0) {
        return {std::move(tour), std::numeric_limits<double>::max()};
    }

    greedy_insert(prob, tour, in_tour, remaining_cap, order, edge_active);

    // Ensure enough customers for a valid solution.
    // Tour needs >= 4 elements [depot, a, b, depot] (binary edges: can't use same edge twice).
    // Path needs >= 2 elements [source, target] (always valid, edge source->target).
    int32_t min_valid_size = is_tour ? 4 : 2;
    int32_t min_construction_size = is_tour ? 3 : 2;

    while (static_cast<int32_t>(tour.size()) <= min_construction_size && !customers.empty()) {
        double bd = std::numeric_limits<double>::max();
        int32_t bn = -1;
        for (int32_t c : customers) {
            if (in_tour[c]) continue;
            if (prob.demand(c) > remaining_cap) continue;
            double best_delta = std::numeric_limits<double>::max();
            for (size_t pos = 0; pos + 1 < tour.size(); ++pos) {
                double delta = edge_cost(prob, tour[pos], c, edge_active)
                             + edge_cost(prob, c, tour[pos + 1], edge_active)
                             - edge_cost(prob, tour[pos], tour[pos + 1], edge_active)
                             - prob.profit(c);
                if (delta < best_delta) best_delta = delta;
            }
            if (best_delta < bd) { bd = best_delta; bn = c; }
        }
        if (bn < 0) break;
        size_t best_pos = 1;
        double best_delta = std::numeric_limits<double>::max();
        for (size_t pos = 0; pos + 1 < tour.size(); ++pos) {
            double delta = edge_cost(prob, tour[pos], bn, edge_active)
                         + edge_cost(prob, bn, tour[pos + 1], edge_active)
                         - edge_cost(prob, tour[pos], tour[pos + 1], edge_active);
            if (delta < best_delta) { best_delta = delta; best_pos = pos + 1; }
        }
        tour.insert(tour.begin() + static_cast<ptrdiff_t>(best_pos), bn);
        in_tour[bn] = true;
        remaining_cap -= prob.demand(bn);
    }

    if (tour.size() > 2) {
        local_search(prob, tour, in_tour, remaining_cap, edge_active);
    }

    double obj = (static_cast<int32_t>(tour.size()) >= min_valid_size)
                     ? tour_objective(prob, tour)
                     : std::numeric_limits<double>::max();
    return {std::move(tour), obj};
}

/// Convert a tour (node sequence) to a MIP solution vector.
inline std::vector<double> tour_to_solution(const Problem& prob,
                                            const std::vector<int32_t>& tour) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    int32_t min_valid_size = prob.is_tour() ? 4 : 2;

    std::vector<double> sol(m + n, 0.0);
    sol[m + source] = 1.0;
    sol[m + target] = 1.0;

    if (static_cast<int32_t>(tour.size()) >= min_valid_size) {
        for (size_t i = 1; i + 1 < tour.size(); ++i)
            sol[m + tour[i]] = 1.0;
        for (size_t i = 0; i + 1 < tour.size(); ++i) {
            int32_t e = find_edge(g, tour[i], tour[i + 1]);
            if (e >= 0) sol[e] = 1.0;
        }
    }
    return sol;
}

}  // namespace detail

struct HeuristicResult {
    std::vector<double> col_values;  // size (num_edges + num_nodes)
    double objective;                // minimization objective (cost - profit)
};

// ─── Graph reduction strategies ───────────────────────────────────────────

/// Strategy A: LP-value threshold.
/// Include edges with x_e > threshold + edges between pairs of nodes with y_i > node_threshold.
inline ReducedGraph reduce_lp_threshold(const Problem& prob,
                                        std::span<const double> x_lp,
                                        std::span<const double> y_lp,
                                        double edge_threshold = 0.1,
                                        double node_threshold = 0.5) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    ReducedGraph rg;
    rg.edge_active.assign(m, false);
    rg.node_active.assign(n, false);

    // Always activate source and target
    rg.node_active[prob.source()] = true;
    rg.node_active[prob.target()] = true;

    // Activate nodes with high y-value
    for (int32_t i = 0; i < n; ++i) {
        if (y_lp[i] > node_threshold) rg.node_active[i] = true;
    }

    // Activate edges above threshold or incident to active nodes
    for (auto e : g.edges()) {
        int32_t u = g.edge_source(e);
        int32_t v = g.edge_target(e);
        if (x_lp[e] > edge_threshold
            || (rg.node_active[u] && rg.node_active[v])) {
            rg.edge_active[e] = true;
            rg.node_active[u] = true;
            rg.node_active[v] = true;
        }
    }

    return rg;
}

/// Strategy B: RINS-style (when incumbent exists).
/// Include edges from incumbent + fractional LP edges.
/// No neighborhood expansion (complete graphs make that useless).
inline ReducedGraph reduce_rins(const Problem& prob,
                                std::span<const double> x_lp,
                                std::span<const double> y_lp,
                                std::span<const double> incumbent,
                                double lp_threshold = 0.1) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    // Fall back to Strategy A if no incumbent
    if (incumbent.empty() || static_cast<int>(incumbent.size()) < m + n) {
        return reduce_lp_threshold(prob, x_lp, y_lp);
    }

    ReducedGraph rg;
    rg.edge_active.assign(m, false);
    rg.node_active.assign(n, false);

    rg.node_active[prob.source()] = true;
    rg.node_active[prob.target()] = true;

    // Include edges in incumbent OR fractional in LP
    for (auto e : g.edges()) {
        if (incumbent[e] > 0.5 || x_lp[e] > lp_threshold) {
            rg.edge_active[e] = true;
            rg.node_active[g.edge_source(e)] = true;
            rg.node_active[g.edge_target(e)] = true;
        }
    }

    // Activate nodes from incumbent or LP
    for (int32_t i = 0; i < n; ++i) {
        if (incumbent[m + i] > 0.5 || y_lp[i] > 0.5)
            rg.node_active[i] = true;
    }

    return rg;
}

/// Strategy C: Neighborhood expansion.
/// Seed edges with x_e > seed_threshold, then expand 1-hop.
inline ReducedGraph reduce_neighborhood(const Problem& prob,
                                        std::span<const double> x_lp,
                                        std::span<const double> y_lp,
                                        double seed_threshold = 0.3) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    ReducedGraph rg;
    rg.edge_active.assign(m, false);
    rg.node_active.assign(n, false);

    rg.node_active[prob.source()] = true;
    rg.node_active[prob.target()] = true;

    // Seed nodes: endpoints of edges with x_e > threshold
    for (auto e : g.edges()) {
        if (x_lp[e] > seed_threshold) {
            rg.node_active[g.edge_source(e)] = true;
            rg.node_active[g.edge_target(e)] = true;
        }
    }

    // Activate seed edges + all edges between active node pairs
    for (auto e : g.edges()) {
        int32_t u = g.edge_source(e);
        int32_t v = g.edge_target(e);
        if (x_lp[e] > seed_threshold
            || (rg.node_active[u] && rg.node_active[v])) {
            rg.edge_active[e] = true;
        }
    }

    return rg;
}

// ─── Public heuristic functions ───────────────────────────────────────────

/// Build an initial solution via parallel randomized construction + local search.
/// Runs on the complete graph. Used before MIP solve.
/// num_restarts: total number of restarts (including 3 deterministic orderings).
/// time_budget_ms: when > 0, use opportunistic (non-deterministic) time-based
///     restarts instead of a fixed count — allows more restarts on fast hardware.
/// Returns solution vector + objective value.
inline HeuristicResult build_initial_solution(const Problem& prob,
                                              int num_restarts = 50,
                                              double time_budget_ms = 0.0) {
    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();

    // Build candidate list (exclude source and target)
    std::vector<int32_t> customers;
    for (int32_t i = 0; i < n; ++i) {
        if (i != source && i != target && prob.demand(i) <= Q)
            customers.push_back(i);
    }

    // Pre-build deterministic orderings
    std::vector<std::vector<int32_t>> fixed_orders;

    // Order 1: profit/demand ratio (descending)
    {
        auto order = customers;
        std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            double ra = prob.demand(a) > 0 ? prob.profit(a) / prob.demand(a)
                                           : prob.profit(a);
            double rb = prob.demand(b) > 0 ? prob.profit(b) / prob.demand(b)
                                           : prob.profit(b);
            return ra > rb;
        });
        fixed_orders.push_back(std::move(order));
    }
    // Order 2: profit (descending)
    {
        auto order = customers;
        std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            return prob.profit(a) > prob.profit(b);
        });
        fixed_orders.push_back(std::move(order));
    }
    // Order 3: cheapest from source
    {
        auto order = customers;
        std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            return detail::edge_cost(prob, source, a)
                 < detail::edge_cost(prob, source, b);
        });
        fixed_orders.push_back(std::move(order));
    }

    std::vector<int32_t> best_tour;
    double best_obj = std::numeric_limits<double>::max();

    if (time_budget_ms > 0.0) {
        // --- Opportunistic mode: time-based workers, non-deterministic seeds ---
        auto deadline = std::chrono::steady_clock::now()
                      + std::chrono::microseconds(
                            static_cast<int64_t>(time_budget_ms * 1000));

        std::mutex best_mtx;
        std::atomic<int> next_restart{0};

        auto worker = [&]() {
            std::mt19937 rng(std::random_device{}());
            while (std::chrono::steady_clock::now() < deadline) {
                int r = next_restart.fetch_add(1, std::memory_order_relaxed);
                std::vector<int32_t> order;
                if (r < static_cast<int>(fixed_orders.size())) {
                    order = fixed_orders[r];
                } else {
                    order = customers;
                    std::shuffle(order.begin(), order.end(), rng);
                }
                auto [tour, obj] = detail::single_restart(prob, customers, order);
                if (obj < best_obj) {
                    std::lock_guard lock(best_mtx);
                    if (obj < best_obj) {
                        best_obj = obj;
                        best_tour = std::move(tour);
                    }
                }
            }
        };

        unsigned hw = std::thread::hardware_concurrency();
        unsigned num_workers = std::max(1u, hw);
        {
            tbb::task_group tg;
            for (unsigned i = 0; i < num_workers; ++i)
                tg.run(worker);
            tg.wait();
        }
    } else {
        // --- Deterministic mode: fixed restarts, deterministic seeds ---
        const int num_fixed = static_cast<int>(fixed_orders.size());
        const int num_random = std::max(0, num_restarts - num_fixed);
        fixed_orders.reserve(num_fixed + num_random);
        for (int r = 0; r < num_random; ++r) {
            std::mt19937 rng(static_cast<uint32_t>(r));
            auto order = customers;
            std::shuffle(order.begin(), order.end(), rng);
            fixed_orders.push_back(std::move(order));
        }

        const int total_restarts = static_cast<int>(fixed_orders.size());

        struct RestartResult {
            std::vector<int32_t> tour;
            double obj;
        };
        std::vector<RestartResult> results(total_restarts);

        tbb::parallel_for(0, total_restarts, [&](int i) {
            auto [tour, obj] = detail::single_restart(
                prob, customers, fixed_orders[i]);
            results[i] = {std::move(tour), obj};
        });

        for (int i = 0; i < total_restarts; ++i) {
            if (results[i].obj < best_obj) {
                best_obj = results[i].obj;
                best_tour = std::move(results[i].tour);
            }
        }
    }

    // Convert best tour to MIP solution vector
    auto sol = detail::tour_to_solution(prob, best_tour);
    return {std::move(sol), best_obj};
}

/// Compatibility wrapper: build_warm_start delegates to build_initial_solution.
inline HeuristicResult build_warm_start(const Problem& prob,
                                        int num_restarts = 50,
                                        double time_budget_ms = 0.0) {
    return build_initial_solution(prob, num_restarts, time_budget_ms);
}

/// LP-guided primal heuristic: builds reduced graphs from LP relaxation
/// values and runs combinatorial construction + local search on them.
/// strategy: 0 = all (default), 1 = LP-threshold only, 2 = RINS only, 3 = neighborhood only
/// Returns improved solution, or empty col_values if none found.
inline HeuristicResult lp_guided_heuristic(
    const Problem& prob,
    std::span<const double> x_lp,
    std::span<const double> y_lp,
    std::span<const double> incumbent,
    double time_budget_ms = 20.0,
    int strategy = 0) {

    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();

    auto deadline = std::chrono::steady_clock::now()
                  + std::chrono::microseconds(static_cast<int64_t>(time_budget_ms * 1000));

    // Build reduced graphs for selected strategies
    std::vector<ReducedGraph> graphs;
    if (strategy == 0 || strategy == 1)
        graphs.push_back(reduce_lp_threshold(prob, x_lp, y_lp));
    if (strategy == 0 || strategy == 2)
        graphs.push_back(reduce_rins(prob, x_lp, y_lp, incumbent));
    if (strategy == 0 || strategy == 3)
        graphs.push_back(reduce_neighborhood(prob, x_lp, y_lp));

    std::vector<const ReducedGraph*> strategies;
    for (auto& g : graphs) strategies.push_back(&g);
    const int num_strategies = static_cast<int>(strategies.size());
    if (num_strategies == 0)
        return {{}, std::numeric_limits<double>::max()};

    // Build customer lists per strategy (only active nodes)
    std::vector<std::vector<int32_t>> strategy_customers(num_strategies);
    for (int s = 0; s < num_strategies; ++s) {
        for (int32_t i = 0; i < n; ++i) {
            if (i != source && i != target
                && strategies[s]->node_active[i]
                && prob.demand(i) <= Q) {
                strategy_customers[s].push_back(i);
            }
        }
    }

    // LP-weighted ordering: sort by y_lp descending, breaking ties by profit/demand
    auto make_lp_order = [&](const std::vector<int32_t>& customers) {
        auto order = customers;
        std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            if (std::abs(y_lp[a] - y_lp[b]) > 1e-6) return y_lp[a] > y_lp[b];
            double ra = prob.demand(a) > 0 ? prob.profit(a) / prob.demand(a)
                                           : prob.profit(a);
            double rb = prob.demand(b) > 0 ? prob.profit(b) / prob.demand(b)
                                           : prob.profit(b);
            return ra > rb;
        });
        return order;
    };

    // Shared best solution
    std::mutex best_mtx;
    std::vector<int32_t> best_tour;
    double best_obj = std::numeric_limits<double>::max();
    std::atomic<int> next_restart{0};

    auto worker = [&]() {
        std::mt19937 rng(std::random_device{}());

        while (std::chrono::steady_clock::now() < deadline) {
            int r = next_restart.fetch_add(1, std::memory_order_relaxed);
            int strat = r % num_strategies;

            const auto& customers = strategy_customers[strat];
            const auto& ea = strategies[strat]->edge_active;

            std::vector<int32_t> order;
            if (r < num_strategies) {
                // First round per strategy: LP-weighted ordering
                order = make_lp_order(customers);
            } else {
                // Subsequent rounds: random shuffle
                order = customers;
                std::shuffle(order.begin(), order.end(), rng);
            }

            auto [tour, obj] = detail::single_restart(prob, customers, order, ea);

            if (obj < best_obj) {
                std::lock_guard lock(best_mtx);
                if (obj < best_obj) {
                    best_obj = obj;
                    best_tour = std::move(tour);
                }
            }
        }
    };

    // Launch parallel workers
    unsigned hw = std::thread::hardware_concurrency();
    unsigned num_workers = std::max(1u, hw);

    {
        tbb::task_group tg;
        for (unsigned i = 0; i < num_workers; ++i)
            tg.run(worker);
        tg.wait();
    }

    if (best_tour.empty()) {
        return {{}, std::numeric_limits<double>::max()};
    }

    auto sol = detail::tour_to_solution(prob, best_tour);
    return {std::move(sol), best_obj};
}

}  // namespace rcspp::heuristic
