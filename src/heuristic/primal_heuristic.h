#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <span>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include "core/problem.h"
#include "preprocess/edge_elimination.h"
namespace cptp::heuristic {

/// Subset of the original problem's edges/nodes, built from LP relaxation values.
struct ReducedGraph {
    std::vector<bool> edge_active;  // size m: which edges are in the reduced graph
    std::vector<bool> node_active;  // size n: which nodes are candidates
};

namespace detail {

inline constexpr double kNoSolution = std::numeric_limits<double>::infinity();

struct TourCandidate {
    std::vector<int32_t> tour;
    double objective = kNoSolution;
};

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
    // Minimum size before allowing node drop:
    // tours keep at least [depot, customer, depot] (size 3),
    // paths keep at least [source, target] plus one optional customer (size 3).
    const int32_t min_drop_size = prob.is_tour() ? 4 : 4;

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

inline bool seed_uses_active_edges(const Problem& prob,
                                   std::span<const int32_t> tour,
                                   const std::vector<bool>& edge_active = {}) {
    constexpr double kBig = 1e15;
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        if (edge_cost(prob, tour[i], tour[i + 1], edge_active) > kBig) {
            return false;
        }
    }
    return true;
}

inline bool init_seed_state(const Problem& prob,
                            std::span<const int32_t> tour,
                            std::vector<bool>& in_tour,
                            double& remaining_cap) {
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    in_tour.assign(static_cast<size_t>(n), false);
    for (int32_t v : tour) {
        if (v < 0 || v >= n) return false;
        in_tour[v] = true;
    }

    remaining_cap = prob.capacity() - prob.demand(source);
    if (!prob.is_tour()) remaining_cap -= prob.demand(target);
    if (remaining_cap < -1e-9) return false;

    for (size_t i = 1; i + 1 < tour.size(); ++i) {
        const int32_t v = tour[i];
        if (v != source && v != target) remaining_cap -= prob.demand(v);
    }
    return remaining_cap >= -1e-9;
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
        return {std::move(tour), kNoSolution};
    }

    greedy_insert(prob, tour, in_tour, remaining_cap, order, edge_active);

    // Ensure enough nodes for a valid solution.
    // Tour allows [depot, customer, depot] because depot-incident x_e can be 2.
    // Path needs >= 2 elements [source, target].
    int32_t min_valid_size = is_tour ? 3 : 2;
    int32_t min_construction_size = 2;

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
                     : kNoSolution;
    return {std::move(tour), obj};
}

/// Deterministic fallback constructor used when restart search finds no feasible tour/path.
inline std::vector<int32_t> fallback_feasible(const Problem& prob,
                                              const std::vector<bool>& edge_active = {}) {
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const bool is_tour = prob.is_tour();

    double remaining_cap = prob.capacity() - prob.demand(source);
    if (!is_tour) remaining_cap -= prob.demand(target);
    if (remaining_cap < 0.0) return {};

    constexpr double kBig = 1e15;
    if (is_tour) {
        int32_t best = -1;
        double best_obj = kNoSolution;
        for (int32_t c = 0; c < n; ++c) {
            if (c == source) continue;
            if (prob.demand(c) > remaining_cap) continue;
            const double c1 = edge_cost(prob, source, c, edge_active);
            const double c2 = edge_cost(prob, c, source, edge_active);
            if (c1 > kBig || c2 > kBig) continue;
            const std::vector<int32_t> cand = {source, c, source};
            const double obj = tour_objective(prob, cand);
            if (obj < best_obj) {
                best_obj = obj;
                best = c;
            }
        }
        if (best >= 0) return {source, best, source};
        return {};
    }

    const double direct = edge_cost(prob, source, target, edge_active);
    if (direct <= kBig) return {source, target};

    int32_t best = -1;
    double best_obj = kNoSolution;
    for (int32_t c = 0; c < n; ++c) {
        if (c == source || c == target) continue;
        if (prob.demand(c) > remaining_cap) continue;
        const double c1 = edge_cost(prob, source, c, edge_active);
        const double c2 = edge_cost(prob, c, target, edge_active);
        if (c1 > kBig || c2 > kBig) continue;
        const std::vector<int32_t> cand = {source, c, target};
        const double obj = tour_objective(prob, cand);
        if (obj < best_obj) {
            best_obj = obj;
            best = c;
        }
    }
    if (best >= 0) return {source, best, target};
    return {};
}

/// Enumerate all 1-customer / 2-customer seeds and return sorted candidates.
inline std::vector<TourCandidate>
enumerate_small_seed_candidates(const Problem& prob,
                                const std::vector<int32_t>& customers,
                                const std::vector<bool>& edge_active = {}) {
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const bool is_tour = prob.is_tour();

    double base_remaining = prob.capacity() - prob.demand(source);
    if (!is_tour) base_remaining -= prob.demand(target);
    if (base_remaining < 0.0) return {};

    std::vector<TourCandidate> candidates;
    candidates.reserve(static_cast<size_t>(customers.size()) * customers.size());

    auto try_seed = [&](std::vector<int32_t> seed, double used_demand) {
        if (used_demand > base_remaining + 1e-9) return;
        if (!seed_uses_active_edges(prob, seed, edge_active)) return;
        candidates.push_back(TourCandidate{
            .tour = std::move(seed),
            .objective = kNoSolution
        });
    };

    // Path-only zero-customer seed.
    if (!is_tour) {
        try_seed({source, target}, 0.0);
    }

    // All single-customer seeds.
    for (int32_t c : customers) {
        double used = prob.demand(c);
        if (is_tour) {
            try_seed({source, c, source}, used);
        } else {
            try_seed({source, c, target}, used);
        }
    }

    // All two-customer seeds.
    for (size_t i = 0; i < customers.size(); ++i) {
        const int32_t a = customers[i];
        for (size_t j = i + 1; j < customers.size(); ++j) {
            const int32_t b = customers[j];
            const double used = prob.demand(a) + prob.demand(b);
            if (is_tour) {
                try_seed({source, a, b, source}, used);
                try_seed({source, b, a, source}, used);
            } else {
                try_seed({source, a, b, target}, used);
                try_seed({source, b, a, target}, used);
            }
        }
    }

    for (auto& cand : candidates) {
        cand.objective = tour_objective(prob, cand.tour);
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const TourCandidate& a, const TourCandidate& b) {
                         if (std::abs(a.objective - b.objective) > 1e-12) {
                             return a.objective < b.objective;
                         }
                         return a.tour < b.tour;
                     });
    return candidates;
}

/// Convert a tour (node sequence) to a MIP solution vector.
inline std::vector<double> tour_to_solution(const Problem& prob,
                                            const std::vector<int32_t>& tour) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    int32_t min_valid_size = prob.is_tour() ? 3 : 2;

    std::vector<double> sol(m + n, 0.0);
    sol[m + source] = 1.0;
    sol[m + target] = 1.0;

    if (static_cast<int32_t>(tour.size()) >= min_valid_size) {
        for (size_t i = 1; i + 1 < tour.size(); ++i)
            sol[m + tour[i]] = 1.0;
        for (size_t i = 0; i + 1 < tour.size(); ++i) {
            int32_t e = find_edge(g, tour[i], tour[i + 1]);
            if (e >= 0) sol[e] += 1.0;
        }
    }
    return sol;
}

}  // namespace detail

struct HeuristicResult {
    std::vector<double> col_values;  // size (num_edges + num_nodes)
    double objective;                // minimization objective (cost - profit)
};

struct ConstructionPool {
    std::vector<detail::TourCandidate> candidates;  // sorted by objective
    std::vector<bool> edge_active;                  // optional reduced graph
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

inline std::vector<bool> build_edge_activity_from_bounds(
    const Problem& prob,
    std::span<const double> fwd_bounds,
    std::span<const double> bwd_bounds,
    double correction,
    double upper_bound) {
    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();
    if (!fwd_bounds.empty() && !bwd_bounds.empty()
        && fwd_bounds.size() == static_cast<size_t>(n)
        && bwd_bounds.size() == static_cast<size_t>(n)
        && std::isfinite(upper_bound)) {
        std::vector<double> fwd_vec(fwd_bounds.begin(), fwd_bounds.end());
        std::vector<double> bwd_vec(bwd_bounds.begin(), bwd_bounds.end());
        auto eliminated = preprocess::edge_elimination(
            prob, fwd_vec, bwd_vec, upper_bound, correction);
        std::vector<bool> edge_active(static_cast<size_t>(m), true);
        int32_t inactive = 0;
        for (int32_t e = 0; e < m; ++e) {
            if (eliminated[e]) {
                edge_active[e] = false;
                inactive++;
            }
        }
        if (inactive > 0) return edge_active;
    }
    return {};
}

inline ConstructionPool build_construction_pool(
    const Problem& prob,
    int max_candidates = 0,
    std::span<const double> fwd_bounds = {},
    std::span<const double> bwd_bounds = {},
    double correction = 0.0,
    double upper_bound = std::numeric_limits<double>::infinity()) {
    ConstructionPool pool;
    pool.edge_active = build_edge_activity_from_bounds(
        prob, fwd_bounds, bwd_bounds, correction, upper_bound);

    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();
    std::vector<int32_t> customers;
    for (int32_t i = 0; i < n; ++i) {
        if (i != source && i != target && prob.demand(i) <= Q) {
            customers.push_back(i);
        }
    }

    pool.candidates = detail::enumerate_small_seed_candidates(
        prob, customers, pool.edge_active);
    if (pool.candidates.empty()) {
        auto fallback = detail::fallback_feasible(prob, pool.edge_active);
        if (!fallback.empty()) {
            const double fallback_obj = detail::tour_objective(prob, fallback);
            pool.candidates.push_back(detail::TourCandidate{
                .tour = std::move(fallback),
                .objective = fallback_obj
            });
        }
    }
    if (max_candidates > 0
        && static_cast<int32_t>(pool.candidates.size()) > max_candidates) {
        pool.candidates.resize(static_cast<size_t>(max_candidates));
    }
    return pool;
}

inline HeuristicResult best_construction_solution(const Problem& prob,
                                                  const ConstructionPool& pool) {
    if (pool.candidates.empty()) {
        auto fallback = detail::fallback_feasible(prob, pool.edge_active);
        if (fallback.empty()) return {{}, detail::kNoSolution};
        const double obj = detail::tour_objective(prob, fallback);
        return {detail::tour_to_solution(prob, fallback), obj};
    }

    const auto& best = pool.candidates.front();
    return {detail::tour_to_solution(prob, best.tour), best.objective};
}

inline HeuristicResult run_local_search_from_pool(
    const Problem& prob,
    const ConstructionPool& pool,
    int num_starts,
    int max_ls_iter = 200,
    int32_t max_workers = 0,
    const std::vector<bool>& edge_active_override = {}) {
    if (pool.candidates.empty()) {
        return {{}, detail::kNoSolution};
    }

    const std::vector<bool>& edge_active =
        edge_active_override.empty() ? pool.edge_active : edge_active_override;
    int starts = std::clamp(num_starts, 1, static_cast<int>(pool.candidates.size()));
    if (starts <= 0) {
        return {{}, detail::kNoSolution};
    }

    std::vector<detail::TourCandidate> improved(static_cast<size_t>(starts));
    auto run_ls = [&]() {
        tbb::parallel_for(0, starts, [&](int i) {
            const auto& cand = pool.candidates[static_cast<size_t>(i)];
            std::vector<int32_t> tour = cand.tour;
            if (!detail::seed_uses_active_edges(prob, tour, edge_active)) {
                improved[static_cast<size_t>(i)] = detail::TourCandidate{};
                return;
            }

            std::vector<bool> in_tour;
            double remaining_cap = 0.0;
            if (!detail::init_seed_state(prob, tour, in_tour, remaining_cap)) {
                improved[static_cast<size_t>(i)] = detail::TourCandidate{};
                return;
            }

            detail::local_search(prob, tour, in_tour, remaining_cap,
                                 edge_active, max_ls_iter);
            const double obj = detail::tour_objective(prob, tour);
            improved[static_cast<size_t>(i)] = detail::TourCandidate{
                .tour = std::move(tour),
                .objective = obj
            };
        });
    };
    if (max_workers > 0) {
        tbb::task_arena arena(max_workers);
        arena.execute(run_ls);
    } else {
        run_ls();
    }

    detail::TourCandidate best;
    for (const auto& cand : improved) {
        if (cand.tour.empty()) continue;
        if (!std::isfinite(best.objective)
            || cand.objective + 1e-12 < best.objective
            || (std::abs(cand.objective - best.objective) <= 1e-12
                && cand.tour < best.tour)) {
            best = cand;
        }
    }
    if (best.tour.empty()) {
        return {{}, detail::kNoSolution};
    }
    return {detail::tour_to_solution(prob, best.tour), best.objective};
}

/// Build an initial solution via deterministic seed construction + parallel local
/// search over top seeds.
inline HeuristicResult build_initial_solution(const Problem& prob,
                                              int num_restarts = 50,
                                              double time_budget_ms = 0.0,
                                              std::span<const double> fwd_bounds = {},
                                              std::span<const double> bwd_bounds = {},
                                              double correction = 0.0,
                                              double upper_bound = std::numeric_limits<double>::infinity(),
                                              int32_t max_workers = 0) {
    (void)time_budget_ms;  // deterministic-only workflow

    const int32_t candidate_cap = std::max<int32_t>(
        1, std::max<int32_t>(num_restarts, (max_workers > 0 ? max_workers : 1)));
    auto pool = build_construction_pool(
        prob, candidate_cap, fwd_bounds, bwd_bounds, correction, upper_bound);
    auto construction = best_construction_solution(prob, pool);

    const int32_t ls_starts = std::min<int32_t>(
        candidate_cap,
        std::max<int32_t>(1, max_workers > 0 ? max_workers : num_restarts));
    auto improved = run_local_search_from_pool(
        prob, pool, ls_starts, 200, max_workers);

    if (improved.col_values.empty()
        || improved.objective > construction.objective + 1e-12) {
        return construction;
    }
    return improved;
}

/// Compatibility wrapper: build_warm_start delegates to build_initial_solution.
inline HeuristicResult build_warm_start(const Problem& prob,
                                        int num_restarts = 50,
                                        double time_budget_ms = 0.0,
                                        std::span<const double> fwd_bounds = {},
                                        std::span<const double> bwd_bounds = {},
                                        double correction = 0.0,
                                        double upper_bound = std::numeric_limits<double>::infinity(),
                                        int32_t max_workers = 0) {
    return build_initial_solution(
        prob, num_restarts, time_budget_ms,
        fwd_bounds, bwd_bounds, correction, upper_bound, max_workers);
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
    int strategy = 0,
    int max_restarts = 0,
    uint32_t restart_seed = 0,
    double edge_threshold = 0.1,
    double node_threshold = 0.5,
    double lp_threshold = 0.1,
    double seed_threshold = 0.3) {

    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();
    (void)time_budget_ms;  // deterministic-only workflow

    // Build reduced graphs for selected strategies
    std::vector<ReducedGraph> graphs;
    if (strategy == 0 || strategy == 1)
        graphs.push_back(reduce_lp_threshold(prob, x_lp, y_lp, edge_threshold, node_threshold));
    if (strategy == 0 || strategy == 2)
        graphs.push_back(reduce_rins(prob, x_lp, y_lp, incumbent, lp_threshold));
    if (strategy == 0 || strategy == 3)
        graphs.push_back(reduce_neighborhood(prob, x_lp, y_lp, seed_threshold));

    std::vector<const ReducedGraph*> strategies;
    for (auto& g : graphs) strategies.push_back(&g);
    const int num_strategies = static_cast<int>(strategies.size());
    if (num_strategies == 0)
        return {{}, detail::kNoSolution};

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

    const int restart_target =
        (max_restarts > 0) ? max_restarts : std::max(1, 8 * num_strategies);
    int total_restarts = restart_target;
    if (total_restarts <= 0) {
        return {{}, detail::kNoSolution};
    }

    struct RestartResult {
        std::vector<int32_t> tour;
        double obj = detail::kNoSolution;
    };
    std::vector<RestartResult> results(static_cast<size_t>(total_restarts));

    tbb::parallel_for(0, total_restarts, [&](int r) {
        const int strat = r % num_strategies;
        const auto& customers = strategy_customers[strat];
        const auto& ea = strategies[strat]->edge_active;

        std::vector<int32_t> order;
        if (r < num_strategies) {
            order = make_lp_order(customers);
        } else {
            order = customers;
            std::mt19937 rng(restart_seed + static_cast<uint32_t>(r));
            std::shuffle(order.begin(), order.end(), rng);
        }

        auto [tour, obj] = detail::single_restart(prob, customers, order, ea);
        results[static_cast<size_t>(r)] = {std::move(tour), obj};
    });

    double best_obj = detail::kNoSolution;
    std::vector<int32_t> best_tour;
    for (int r = 0; r < total_restarts; ++r) {
        const auto& rr = results[static_cast<size_t>(r)];
        if (rr.obj < best_obj) {
            best_obj = rr.obj;
            best_tour = rr.tour;
        }
    }
    if (best_tour.empty()) {
        return {{}, detail::kNoSolution};
    }
    auto sol = detail::tour_to_solution(prob, best_tour);
    return {std::move(sol), best_obj};
}

}  // namespace cptp::heuristic
