#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include "core/problem.h"
#include "preprocess/edge_elimination.h"
#include "sep/separation_oracle.h"

namespace cptp::heuristic {

/// Subset of the original problem's edges/nodes, built from LP relaxation values.
struct ReducedGraph {
    std::vector<bool> edge_active;  // size m: which edges are in the reduced graph
    std::vector<bool> node_active;  // size n: which nodes are candidates
};

inline constexpr int32_t kDefaultLsMaxIterPerStart = 1000;

namespace detail {

inline constexpr double kNoSolution = std::numeric_limits<double>::infinity();
inline constexpr double kLsEpsImprove = 1e-9;
inline constexpr double kLsEpsTie = 1e-12;
inline constexpr double kBigEdgeCost = 1e15;
inline constexpr int32_t kIlsRounds = 8;
inline constexpr int32_t kKickTrials = 16;

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
            if (cost_pc > kBigEdgeCost || cost_cn > kBigEdgeCost) continue;

            double cost_pn = edge_cost(prob, prev, next, edge_active);
            if (cost_pn > kBigEdgeCost) cost_pn = 0.0;

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
inline int local_search(const Problem& prob,
                        std::vector<int32_t>& tour,
                        std::vector<bool>& in_tour,
                        double& remaining_cap,
                        const std::vector<bool>& edge_active = {},
                        int max_iter_per_start = kDefaultLsMaxIterPerStart) {
    const int32_t n = prob.num_nodes();
    // Minimum size before allowing node drop:
    // tours keep at least [depot, customer, depot] (size 3),
    // paths keep at least [source, target] plus one optional customer (size 3).
    const int32_t min_drop_size = prob.is_tour() ? 4 : 4;

    int iterations = 0;
    for (int iter = 0; iter < max_iter_per_start; ++iter) {
        iterations++;
        bool improved = false;
        int32_t len = static_cast<int32_t>(tour.size());

        // --- 2-opt ---
        for (int32_t i = 0; i < len - 2 && !improved; ++i) {
            for (int32_t j = i + 2; j < len - 1; ++j) {
                double old_c = edge_cost(prob, tour[i], tour[i + 1], edge_active)
                             + edge_cost(prob, tour[j], tour[j + 1], edge_active);
                double new_c = edge_cost(prob, tour[i], tour[j], edge_active)
                             + edge_cost(prob, tour[i + 1], tour[j + 1], edge_active);
                if (new_c < old_c - kLsEpsImprove) {
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
                    if (insert_add - remove_save < -kLsEpsImprove) {
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

        // --- Swap (1,1): exchange one visited customer with one unvisited customer ---
        for (int32_t i = 1; i < len - 1 && !improved; ++i) {
            const int32_t out = tour[i];
            if (out == prob.source() || out == prob.target()) continue;

            for (int32_t in = 0; in < n && !improved; ++in) {
                if (in_tour[in]) continue;
                if (in == prob.source() || in == prob.target()) continue;

                const double new_remaining =
                    remaining_cap + prob.demand(out) - prob.demand(in);
                if (new_remaining < -kLsEpsImprove) continue;

                const double old_prev = edge_cost(prob, tour[i - 1], out, edge_active);
                const double old_next = edge_cost(prob, out, tour[i + 1], edge_active);
                const double new_prev = edge_cost(prob, tour[i - 1], in, edge_active);
                const double new_next = edge_cost(prob, in, tour[i + 1], edge_active);
                if (old_prev > kBigEdgeCost || old_next > kBigEdgeCost
                    || new_prev > kBigEdgeCost || new_next > kBigEdgeCost) {
                    continue;
                }

                const double delta =
                    (new_prev + new_next - old_prev - old_next)
                    - (prob.profit(in) - prob.profit(out));
                if (delta < -kLsEpsImprove) {
                    tour[i] = in;
                    in_tour[out] = false;
                    in_tour[in] = true;
                    remaining_cap = new_remaining;
                    improved = true;
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
            if (save - prob.profit(c) > kLsEpsImprove) {
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

            double best_delta = kLsEpsImprove;
            size_t best_pos = 0;

            for (size_t pos = 1; pos < tour.size(); ++pos) {
                double cost_pc = edge_cost(prob, tour[pos - 1], c, edge_active);
                double cost_cn = edge_cost(prob, c, tour[pos], edge_active);
                if (cost_pc > kBigEdgeCost || cost_cn > kBigEdgeCost) continue;
                double cost_pn = edge_cost(prob, tour[pos - 1], tour[pos], edge_active);
                if (cost_pn > kBigEdgeCost) cost_pn = 0.0;

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
    return iterations;
}

inline bool seed_uses_active_edges(const Problem& prob,
                                   std::span<const int32_t> tour,
                                   const std::vector<bool>& edge_active = {}) {
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        if (edge_cost(prob, tour[i], tour[i + 1], edge_active) > kBigEdgeCost) {
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
    if (remaining_cap < -kLsEpsImprove) return false;

    for (size_t i = 1; i + 1 < tour.size(); ++i) {
        const int32_t v = tour[i];
        if (v != source && v != target) remaining_cap -= prob.demand(v);
    }
    return remaining_cap >= -kLsEpsImprove;
}

/// SplitMix64 for deterministic index mixing (no RNG state, no time seeding).
inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

inline uint64_t deterministic_kick_key(int32_t start_index,
                                       int32_t round_index,
                                       int32_t trial,
                                       uint64_t salt) {
    uint64_t key = static_cast<uint64_t>(start_index + 1);
    key ^= static_cast<uint64_t>(round_index + 1) * 0x9e3779b97f4a7c15ull;
    key ^= static_cast<uint64_t>(trial + 1) * 0xbf58476d1ce4e5b9ull;
    key ^= salt;
    return mix64(key);
}

inline bool perturb_swap_kick(const Problem& prob,
                              std::vector<int32_t>& tour,
                              std::vector<bool>& in_tour,
                              double& remaining_cap,
                              const std::vector<bool>& edge_active,
                              int32_t start_index,
                              int32_t round_index) {
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const int32_t len = static_cast<int32_t>(tour.size());

    std::vector<int32_t> in_positions;
    in_positions.reserve(tour.size());
    for (int32_t pos = 1; pos + 1 < len; ++pos) {
        const int32_t v = tour[pos];
        if (v == source || v == target) continue;
        in_positions.push_back(pos);
    }
    if (in_positions.empty()) return false;

    std::vector<int32_t> out_nodes;
    out_nodes.reserve(static_cast<size_t>(n));
    for (int32_t c = 0; c < n; ++c) {
        if (in_tour[c]) continue;
        if (c == source || c == target) continue;
        out_nodes.push_back(c);
    }
    if (out_nodes.empty()) return false;

    for (int32_t trial = 0; trial < kKickTrials; ++trial) {
        const uint64_t key = deterministic_kick_key(
            start_index, round_index, trial, 0xd1b54a32d192ed03ull);
        const int32_t pos = in_positions[static_cast<size_t>(
            key % static_cast<uint64_t>(in_positions.size()))];
        const int32_t in = out_nodes[static_cast<size_t>(
            (key >> 32) % static_cast<uint64_t>(out_nodes.size()))];
        const int32_t out = tour[pos];

        const double new_remaining =
            remaining_cap + prob.demand(out) - prob.demand(in);
        if (new_remaining < -kLsEpsImprove) continue;

        const double c1 = edge_cost(prob, tour[pos - 1], in, edge_active);
        const double c2 = edge_cost(prob, in, tour[pos + 1], edge_active);
        if (c1 > kBigEdgeCost || c2 > kBigEdgeCost) continue;

        tour[pos] = in;
        in_tour[out] = false;
        in_tour[in] = true;
        remaining_cap = new_remaining;
        return true;
    }
    return false;
}

inline bool perturb_relocate_kick(const Problem& prob,
                                  std::vector<int32_t>& tour,
                                  const std::vector<bool>& edge_active,
                                  int32_t start_index,
                                  int32_t round_index) {
    const int32_t len = static_cast<int32_t>(tour.size());
    if (len <= 3) return false;
    const int32_t interior_count = len - 2;

    for (int32_t trial = 0; trial < kKickTrials; ++trial) {
        const uint64_t key = deterministic_kick_key(
            start_index, round_index, trial, 0x94d049bb133111ebull);
        const int32_t from = 1 + static_cast<int32_t>(
            key % static_cast<uint64_t>(interior_count));
        const int32_t to = 1 + static_cast<int32_t>(
            (key >> 32) % static_cast<uint64_t>(len - 1));
        if (to == from || to == from + 1) continue;

        std::vector<int32_t> moved = tour;
        const int32_t node = moved[from];
        moved.erase(moved.begin() + from);
        const int32_t ins = (to > from) ? to - 1 : to;
        moved.insert(moved.begin() + ins, node);
        if (!seed_uses_active_edges(prob, moved, edge_active)) continue;

        tour = std::move(moved);
        return true;
    }
    return false;
}

inline bool perturb_deterministic(const Problem& prob,
                                  std::vector<int32_t>& tour,
                                  std::vector<bool>& in_tour,
                                  double& remaining_cap,
                                  const std::vector<bool>& edge_active,
                                  int32_t start_index,
                                  int32_t round_index) {
    if (perturb_swap_kick(prob, tour, in_tour, remaining_cap, edge_active,
                          start_index, round_index)) {
        return true;
    }
    return perturb_relocate_kick(
        prob, tour, edge_active, start_index, round_index);
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
inline void validate_tour_or_throw(const Problem& prob,
                                   std::span<const int32_t> tour) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const int32_t min_valid_size = prob.is_tour() ? 3 : 2;

    if (static_cast<int32_t>(tour.size()) < min_valid_size) {
        throw std::runtime_error("Heuristic produced invalid route: too short");
    }
    if (tour.front() != source || tour.back() != target) {
        throw std::runtime_error("Heuristic produced invalid route: wrong endpoints");
    }

    std::vector<bool> seen(static_cast<size_t>(n), false);
    for (size_t i = 0; i < tour.size(); ++i) {
        const int32_t v = tour[i];
        if (v < 0 || v >= n) {
            throw std::runtime_error("Heuristic produced invalid route: node id out of bounds");
        }
        if (i > 0 && i + 1 < tour.size()) {
            if (v == source || v == target) {
                throw std::runtime_error(
                    "Heuristic produced invalid route: endpoint inside route interior");
            }
            if (seen[static_cast<size_t>(v)]) {
                throw std::runtime_error("Heuristic produced invalid route: repeated customer");
            }
            seen[static_cast<size_t>(v)] = true;
        }
    }

    std::vector<double> x_values(static_cast<size_t>(m), 0.0);
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        const int32_t e = find_edge(g, tour[i], tour[i + 1]);
        if (e < 0) {
            throw std::runtime_error(
                "Heuristic produced invalid route: missing edge in graph");
        }
        x_values[static_cast<size_t>(e)] += 1.0;
    }

    double used = prob.demand(source);
    if (!prob.is_tour()) used += prob.demand(target);
    std::vector<double> y_values(static_cast<size_t>(n), 0.0);
    y_values[static_cast<size_t>(source)] = 1.0;
    y_values[static_cast<size_t>(target)] = 1.0;
    for (int32_t v = 0; v < n; ++v) {
        if (seen[static_cast<size_t>(v)]) {
            used += prob.demand(v);
            y_values[static_cast<size_t>(v)] = 1.0;
        }
    }
    if (used > prob.capacity() + kLsEpsImprove) {
        throw std::runtime_error("Heuristic produced invalid route: capacity violation");
    }

    // Reuse the solver's incumbent SEC validator for consistency.
    sep::SeparationOracle oracle(prob);
    if (!oracle.is_feasible(x_values, y_values, 0, m)) {
        throw std::runtime_error(
            "Heuristic produced invalid route: SEC feasibility violation");
    }
}

inline std::vector<double> tour_to_solution(const Problem& prob,
                                            const std::vector<int32_t>& tour) {
    const auto& g = prob.graph();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    int32_t min_valid_size = prob.is_tour() ? 3 : 2;

    validate_tour_or_throw(prob, tour);

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

struct WarmStartProgressSnapshot {
    int32_t starts_total = 0;
    int32_t starts_done = 0;
    int64_t ls_iterations_total = 0;
    int64_t ub_improvements = 0;
    double best_objective = detail::kNoSolution;
    bool final = false;
};

struct WarmStartProgressOptions {
    int64_t report_every_starts = 8;
    bool report_on_ub_improvement = true;
};

using WarmStartProgressCallback =
    std::function<void(const WarmStartProgressSnapshot&)>;

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
    int ls_max_iter_per_start = kDefaultLsMaxIterPerStart,
    int32_t max_workers = 0,
    const std::vector<bool>& edge_active_override = {},
    const WarmStartProgressOptions* progress_opts = nullptr,
    const WarmStartProgressCallback& progress_cb = {},
    bool reuse_candidates_for_starts = false) {
    if (pool.candidates.empty()) {
        return {{}, detail::kNoSolution};
    }

    const std::vector<bool>& edge_active =
        edge_active_override.empty() ? pool.edge_active : edge_active_override;
    int starts = reuse_candidates_for_starts
        ? std::max(1, num_starts)
        : std::clamp(num_starts, 1, static_cast<int>(pool.candidates.size()));
    if (starts <= 0) {
        return {{}, detail::kNoSolution};
    }

    std::vector<detail::TourCandidate> improved(static_cast<size_t>(starts));
    std::atomic<int32_t> next_start{0};
    std::atomic<int32_t> starts_done{0};
    std::atomic<int64_t> ls_iterations_total{0};
    std::atomic<int64_t> ub_improvements{0};
    const bool progress_enabled = static_cast<bool>(progress_cb) && progress_opts != nullptr;
    const int64_t report_every = progress_enabled
        ? std::max<int64_t>(1, progress_opts->report_every_starts)
        : 0;
    std::atomic<int64_t> next_report_milestone{
        (progress_enabled && report_every > 0)
            ? report_every
            : std::numeric_limits<int64_t>::max()
    };
    std::mutex progress_best_mu;
    double progress_best_objective = detail::kNoSolution;
    std::mutex progress_emit_mu;
    int32_t last_emitted_starts = -1;

    auto build_progress_snapshot = [&](bool final) -> WarmStartProgressSnapshot {
        WarmStartProgressSnapshot snap;
        snap.starts_total = starts;
        snap.starts_done = starts_done.load(std::memory_order_relaxed);
        snap.ls_iterations_total = ls_iterations_total.load(std::memory_order_relaxed);
        snap.ub_improvements = ub_improvements.load(std::memory_order_relaxed);
        {
            std::lock_guard<std::mutex> lock(progress_best_mu);
            snap.best_objective = progress_best_objective;
        }
        snap.final = final;
        return snap;
    };

    auto emit_progress = [&](bool final) {
        if (!progress_enabled) return;
        const auto snap = build_progress_snapshot(final);
        std::lock_guard<std::mutex> lock(progress_emit_mu);
        if (!final && snap.starts_done <= last_emitted_starts) return;
        last_emitted_starts = snap.starts_done;
        progress_cb(snap);
    };

    auto finalize_start = [&](bool ub_improved, int64_t ls_iters_for_start) {
        ls_iterations_total.fetch_add(ls_iters_for_start, std::memory_order_relaxed);
        if (ub_improved) {
            ub_improvements.fetch_add(1, std::memory_order_relaxed);
        }

        const int32_t done_after =
            starts_done.fetch_add(1, std::memory_order_relaxed) + 1;
        if (done_after > starts) {
            throw std::runtime_error(
                "warm-start local search exceeded requested start count");
        }
        bool periodic_due = false;
        if (progress_enabled && report_every > 0) {
            while (true) {
                int64_t milestone =
                    next_report_milestone.load(std::memory_order_relaxed);
                if (done_after < milestone) break;
                if (next_report_milestone.compare_exchange_weak(
                        milestone, milestone + report_every,
                        std::memory_order_relaxed,
                        std::memory_order_relaxed)) {
                    periodic_due = true;
                    break;
                }
            }
        }
        if ((progress_enabled && progress_opts->report_on_ub_improvement && ub_improved)
            || periodic_due) {
            emit_progress(false);
        }
    };

    auto run_ls = [&]() {
        tbb::parallel_for(0, max_workers > 0 ? max_workers : starts, [&](int) {
            while (true) {
                const int32_t i = next_start.fetch_add(1, std::memory_order_relaxed);
                if (i >= starts) break;

                bool ub_improved = false;
                int64_t ls_iters_for_start = 0;
                const int32_t seed_index = reuse_candidates_for_starts
                    ? (i % static_cast<int32_t>(pool.candidates.size()))
                    : i;
                const auto& cand = pool.candidates[static_cast<size_t>(seed_index)];
                std::vector<int32_t> tour = cand.tour;
                if (!detail::seed_uses_active_edges(prob, tour, edge_active)) {
                    improved[static_cast<size_t>(i)] = detail::TourCandidate{};
                    finalize_start(false, ls_iters_for_start);
                    continue;
                }

                std::vector<bool> in_tour;
                double remaining_cap = 0.0;
                if (!detail::init_seed_state(prob, tour, in_tour, remaining_cap)) {
                    improved[static_cast<size_t>(i)] = detail::TourCandidate{};
                    finalize_start(false, ls_iters_for_start);
                    continue;
                }

                // Deterministic contract:
                // - fixed move order in local_search()
                // - fixed number of ILS rounds
                // - deterministic kick candidate order derived from (start, round)
                ls_iters_for_start += detail::local_search(
                    prob, tour, in_tour, remaining_cap, edge_active,
                    ls_max_iter_per_start);
                detail::TourCandidate best_local{
                    .tour = tour,
                    .objective = detail::tour_objective(prob, tour)
                };

                for (int32_t round = 0; round < detail::kIlsRounds; ++round) {
                    std::vector<int32_t> kicked = best_local.tour;
                    std::vector<bool> kicked_in_tour;
                    double kicked_remaining_cap = 0.0;
                    if (!detail::init_seed_state(
                            prob, kicked, kicked_in_tour, kicked_remaining_cap)) {
                        break;
                    }
                    if (!detail::perturb_deterministic(
                            prob, kicked, kicked_in_tour, kicked_remaining_cap,
                            edge_active, i, round)) {
                        continue;
                    }
                    // Recompute state from the perturbed route to keep updates robust.
                    if (!detail::init_seed_state(
                            prob, kicked, kicked_in_tour, kicked_remaining_cap)) {
                        continue;
                    }
                    ls_iters_for_start += detail::local_search(
                        prob, kicked, kicked_in_tour, kicked_remaining_cap,
                        edge_active, ls_max_iter_per_start);
                    const double kicked_obj = detail::tour_objective(prob, kicked);
                    if (kicked_obj + detail::kLsEpsTie < best_local.objective
                        || (std::abs(kicked_obj - best_local.objective)
                                <= detail::kLsEpsTie
                            && kicked < best_local.tour)) {
                        best_local.tour = std::move(kicked);
                        best_local.objective = kicked_obj;
                    }
                }

                improved[static_cast<size_t>(i)] = std::move(best_local);
                if (std::isfinite(improved[static_cast<size_t>(i)].objective)) {
                    std::lock_guard<std::mutex> lock(progress_best_mu);
                    if (!std::isfinite(progress_best_objective)
                        || improved[static_cast<size_t>(i)].objective + detail::kLsEpsTie
                            < progress_best_objective) {
                        progress_best_objective = improved[static_cast<size_t>(i)].objective;
                        ub_improved = true;
                    }
                }
                finalize_start(ub_improved, ls_iters_for_start);
            }
        });
    };
    if (max_workers > 0) {
        tbb::task_arena arena(max_workers);
        arena.execute(run_ls);
    } else {
        run_ls();
    }
    emit_progress(true);

    detail::TourCandidate best;
    for (const auto& cand : improved) {
        if (cand.tour.empty()) continue;
        if (!std::isfinite(best.objective)
            || cand.objective + detail::kLsEpsTie < best.objective
            || (std::abs(cand.objective - best.objective) <= detail::kLsEpsTie
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
        prob, pool, ls_starts, kDefaultLsMaxIterPerStart, max_workers);

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
    double lpg_edge_threshold = 0.1,
    double lpg_node_threshold = 0.5,
    double lpg_lp_threshold = 0.1,
    double lpg_seed_threshold = 0.3) {

    const int32_t n = prob.num_nodes();
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const double Q = prob.capacity();
    (void)time_budget_ms;  // deterministic-only workflow

    // Build reduced graphs for selected strategies
    std::vector<ReducedGraph> graphs;
    if (strategy == 0 || strategy == 1)
        graphs.push_back(reduce_lp_threshold(prob, x_lp, y_lp,
                                             lpg_edge_threshold, lpg_node_threshold));
    if (strategy == 0 || strategy == 2)
        graphs.push_back(reduce_rins(prob, x_lp, y_lp, incumbent,
                                     lpg_lp_threshold));
    if (strategy == 0 || strategy == 3)
        graphs.push_back(reduce_neighborhood(prob, x_lp, y_lp,
                                             lpg_seed_threshold));

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
        if (tour.empty() || !std::isfinite(obj)) {
            results[static_cast<size_t>(r)] = {std::move(tour), obj};
            return;
        }

        detail::TourCandidate best_local{
            .tour = std::move(tour),
            .objective = obj
        };

        const int32_t ils_start = static_cast<int32_t>(
            restart_seed + static_cast<uint32_t>(r));
        for (int32_t round = 0; round < detail::kIlsRounds; ++round) {
            std::vector<int32_t> kicked = best_local.tour;
            std::vector<bool> kicked_in_tour;
            double kicked_remaining_cap = 0.0;
            if (!detail::init_seed_state(
                    prob, kicked, kicked_in_tour, kicked_remaining_cap)) {
                break;
            }
            if (!detail::perturb_deterministic(
                    prob, kicked, kicked_in_tour, kicked_remaining_cap,
                    ea, ils_start, round)) {
                continue;
            }
            // Recompute state from perturbed route before local search.
            if (!detail::init_seed_state(
                    prob, kicked, kicked_in_tour, kicked_remaining_cap)) {
                continue;
            }
            detail::local_search(
                prob, kicked, kicked_in_tour, kicked_remaining_cap, ea,
                kDefaultLsMaxIterPerStart);
            const double kicked_obj = detail::tour_objective(prob, kicked);
            if (kicked_obj + detail::kLsEpsTie < best_local.objective
                || (std::abs(kicked_obj - best_local.objective)
                        <= detail::kLsEpsTie
                    && kicked < best_local.tour)) {
                best_local.tour = std::move(kicked);
                best_local.objective = kicked_obj;
            }
        }

        results[static_cast<size_t>(r)] = {
            std::move(best_local.tour), best_local.objective
        };
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
