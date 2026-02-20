#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <vector>

#include <tbb/task_group.h>

#include "core/problem.h"

namespace cptp::heuristic {

namespace detail {

/// Find edge index between u and v, or -1 if not found.
inline int32_t find_edge(const Graph& g, int32_t u, int32_t v) {
    for (int32_t e : g.incident_edges(u)) {
        if (g.other_endpoint(e, u) == v) return e;
    }
    return -1;
}

inline double edge_cost(const Problem& prob, int32_t u, int32_t v) {
    int32_t e = find_edge(prob.graph(), u, v);
    return (e >= 0) ? prob.edge_cost(e) : std::numeric_limits<double>::max();
}

/// Compute tour objective: sum(edge costs) - sum(profits) - depot_profit.
inline double tour_objective(const Problem& prob,
                             const std::vector<int32_t>& tour) {
    double cost = 0.0;
    for (size_t i = 0; i + 1 < tour.size(); ++i)
        cost += edge_cost(prob, tour[i], tour[i + 1]);
    for (size_t i = 1; i + 1 < tour.size(); ++i)
        cost -= prob.profit(tour[i]);
    cost -= prob.profit(prob.depot());
    return cost;
}

/// Greedy insertion: insert customers in given order at their cheapest position.
inline void greedy_insert(const Problem& prob,
                          std::vector<int32_t>& tour,
                          std::vector<bool>& in_tour,
                          double& remaining_cap,
                          const std::vector<int32_t>& order) {
    for (int32_t c : order) {
        if (in_tour[c]) continue;
        if (prob.demand(c) > remaining_cap) continue;

        double best_delta = std::numeric_limits<double>::max();
        size_t best_pos = 0;

        for (size_t pos = 1; pos < tour.size(); ++pos) {
            int32_t prev = tour[pos - 1];
            int32_t next = tour[pos];

            double cost_pc = edge_cost(prob, prev, c);
            double cost_cn = edge_cost(prob, c, next);
            if (cost_pc > 1e15 || cost_cn > 1e15) continue;

            double cost_pn = edge_cost(prob, prev, next);
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
inline void local_search(const Problem& prob,
                         std::vector<int32_t>& tour,
                         std::vector<bool>& in_tour,
                         double& remaining_cap,
                         int max_iter = 200) {
    const int32_t n = prob.num_nodes();

    for (int iter = 0; iter < max_iter; ++iter) {
        bool improved = false;
        int32_t len = static_cast<int32_t>(tour.size());

        // --- 2-opt ---
        for (int32_t i = 0; i < len - 2 && !improved; ++i) {
            for (int32_t j = i + 2; j < len - 1; ++j) {
                double old_c = edge_cost(prob, tour[i], tour[i + 1])
                             + edge_cost(prob, tour[j], tour[j + 1]);
                double new_c = edge_cost(prob, tour[i], tour[j])
                             + edge_cost(prob, tour[i + 1], tour[j + 1]);
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
                    edge_cost(prob, tour[i - 1], tour[i])
                    + edge_cost(prob, tour[se], tour[se + 1])
                    - edge_cost(prob, tour[i - 1], tour[se + 1]);

                for (int32_t j = 1; j < len - 1 && !improved; ++j) {
                    if (j >= i - 1 && j <= se + 1) continue;
                    double insert_add =
                        edge_cost(prob, tour[j - 1], tour[i])
                        + edge_cost(prob, tour[se], tour[j])
                        - edge_cost(prob, tour[j - 1], tour[j]);
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
        for (int32_t i = 1; i < len - 1 && !improved; ++i) {
            int32_t c = tour[i];
            double save = edge_cost(prob, tour[i - 1], c)
                        + edge_cost(prob, c, tour[i + 1])
                        - edge_cost(prob, tour[i - 1], tour[i + 1]);
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
                double cost_pc = edge_cost(prob, tour[pos - 1], c);
                double cost_cn = edge_cost(prob, c, tour[pos]);
                if (cost_pc > 1e15 || cost_cn > 1e15) continue;
                double cost_pn = edge_cost(prob, tour[pos - 1], tour[pos]);
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
inline std::pair<std::vector<int32_t>, double>
single_restart(const Problem& prob,
               const std::vector<int32_t>& customers,
               const std::vector<int32_t>& order) {
    const int32_t n = prob.num_nodes();
    const int32_t depot = prob.depot();
    const double Q = prob.capacity();

    std::vector<int32_t> tour = {depot, depot};
    std::vector<bool> in_tour(n, false);
    in_tour[depot] = true;
    double remaining_cap = Q;

    greedy_insert(prob, tour, in_tour, remaining_cap, order);

    // Ensure at least one customer
    if (tour.size() <= 2) {
        double bd = std::numeric_limits<double>::max();
        int32_t bn = -1;
        for (int32_t c : customers) {
            double d = edge_cost(prob, depot, c)
                     + edge_cost(prob, c, depot) - prob.profit(c);
            if (d < bd) { bd = d; bn = c; }
        }
        if (bn >= 0) {
            tour.insert(tour.begin() + 1, bn);
            in_tour[bn] = true;
            remaining_cap -= prob.demand(bn);
        }
    }

    if (tour.size() > 2) {
        local_search(prob, tour, in_tour, remaining_cap);
    }

    double obj = (tour.size() > 2) ? tour_objective(prob, tour)
                                   : std::numeric_limits<double>::max();
    return {std::move(tour), obj};
}

}  // namespace detail

/// Build a warm-start solution via parallel randomized construction + local search.
/// time_budget_ms: how long to spend on restarts (default 1000ms).
/// Returns a col_value vector of size (num_edges + num_nodes) for HiGHS setSolution.
inline std::vector<double> build_warm_start(const Problem& prob,
                                            double time_budget_ms = 1000.0) {
    const auto& g = prob.graph();
    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();
    const int32_t depot = prob.depot();
    const double Q = prob.capacity();

    auto deadline = std::chrono::steady_clock::now()
                  + std::chrono::microseconds(static_cast<int64_t>(time_budget_ms * 1000));

    // Build candidate list
    std::vector<int32_t> customers;
    for (int32_t i = 0; i < n; ++i) {
        if (i != depot && prob.demand(i) <= Q)
            customers.push_back(i);
    }

    // Pre-build deterministic orderings
    std::vector<std::vector<int32_t>> fixed_orders;

    // Order 1: profit/demand ratio (descending)
    {
        auto order = customers;
        std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
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
        std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            return prob.profit(a) > prob.profit(b);
        });
        fixed_orders.push_back(std::move(order));
    }
    // Order 3: cheapest from depot
    {
        auto order = customers;
        std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
            return detail::edge_cost(prob, depot, a)
                 < detail::edge_cost(prob, depot, b);
        });
        fixed_orders.push_back(std::move(order));
    }

    // Shared best solution
    std::mutex best_mtx;
    std::vector<int32_t> best_tour;
    double best_obj = std::numeric_limits<double>::max();
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

    // Launch parallel workers
    unsigned hw = std::thread::hardware_concurrency();
    unsigned num_workers = std::max(1u, hw);

    {
        tbb::task_group tg;
        for (unsigned i = 0; i < num_workers; ++i)
            tg.run(worker);
        tg.wait();
    }

    // --- Convert best tour to MIP solution vector ---
    std::vector<double> sol(m + n, 0.0);
    sol[m + depot] = 1.0;

    if (!best_tour.empty() && best_tour.size() > 2) {
        for (size_t i = 1; i + 1 < best_tour.size(); ++i)
            sol[m + best_tour[i]] = 1.0;
        for (size_t i = 0; i + 1 < best_tour.size(); ++i) {
            int32_t e = detail::find_edge(g, best_tour[i], best_tour[i + 1]);
            if (e >= 0) sol[e] = 1.0;
        }
    }

    return sol;
}

}  // namespace cptp::heuristic
