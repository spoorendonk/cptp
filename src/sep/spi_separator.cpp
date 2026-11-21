#include "sep/spi_separator.h"

#include "core/problem.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace cptp::sep {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// ─────────────────────────────────────────────────────────────────────
// Held-Karp DP: exact minimum-cost Hamiltonian path
// ─────────────────────────────────────────────────────────────────────
//
// Given a set S = {nodes[0], ..., nodes[k-1]} of customer nodes, compute
// the minimum over all orderings π of:
//   d[source, π(1)] + d[π(1), π(2)] + ... + d[π(k), target]
//
// The full lower bound on visiting S is:
//   lb = held_karp_result + correction + Σ profit(v) for v in S
//
// dp[mask][last] = min cost path from source through the subset indicated
//                  by `mask`, ending at nodes[last].
//
// Time: O(2^k · k²),  Space: O(2^k · k).

double held_karp(int32_t source, int32_t target, const std::vector<int32_t>& nodes,
                 std::span<const double> ap, int32_t n) {
    const auto k = static_cast<int32_t>(nodes.size());

    if (k == 0) {
        return ap[static_cast<size_t>(source) * n + target];
    }

    auto d = [&](int32_t u, int32_t v) -> double { return ap[static_cast<size_t>(u) * n + v]; };

    if (k == 1) {
        return d(source, nodes[0]) + d(nodes[0], target);
    }

    if (k == 2) {
        double lb1 = d(source, nodes[0]) + d(nodes[0], nodes[1]) + d(nodes[1], target);
        double lb2 = d(source, nodes[1]) + d(nodes[1], nodes[0]) + d(nodes[0], target);
        return std::min(lb1, lb2);
    }

    const int32_t states = 1 << k;
    // Flat array: dp[mask * k + last]
    std::vector<double> dp(static_cast<size_t>(states) * k, kInf);

    // Initialize: single-node subsets
    for (int32_t j = 0; j < k; ++j) {
        dp[static_cast<size_t>(1 << j) * k + j] = d(source, nodes[j]);
    }

    // Fill in order of increasing popcount
    for (int32_t mask = 1; mask < states; ++mask) {
        for (int32_t last = 0; last < k; ++last) {
            if (!(mask & (1 << last))) {
                continue;
            }
            double cur = dp[static_cast<size_t>(mask) * k + last];
            if (cur >= kInf) {
                continue;
            }

            for (int32_t next = 0; next < k; ++next) {
                if (mask & (1 << next)) {
                    continue;
                }
                int32_t new_mask = mask | (1 << next);
                double cost = cur + d(nodes[last], nodes[next]);
                auto& slot = dp[static_cast<size_t>(new_mask) * k + next];
                slot = std::min(cost, slot);
            }
        }
    }

    // Complete: find best way to reach target
    const int32_t full = states - 1;
    double best = kInf;
    for (int32_t last = 0; last < k; ++last) {
        double cost = dp[static_cast<size_t>(full) * k + last] + d(nodes[last], target);
        best = std::min(cost, best);
    }
    return best;
}

/// Compute the lower bound on the objective of any tour/path visiting all
/// nodes in `set`.  Returns -∞ if the set is too large for exact DP.
double compute_set_lb(const Problem& prob, std::span<const double> all_pairs,
                      const std::vector<int32_t>& set, int32_t max_dp_size = 15) {
    const int32_t n = prob.num_nodes();
    const auto k = static_cast<int32_t>(set.size());
    if (k == 0 || k > max_dp_size) {
        return -kInf;
    }

    const int32_t source = prob.source();
    const int32_t target = prob.is_tour() ? source : prob.target();
    const double correction = prob.is_tour() ? prob.profit(source) : 0.0;

    // Junction-node profit correction: each node in the set appears in two
    // consecutive segments, so its profit is subtracted twice.  Add it back.
    double profit_sum = correction;
    for (int32_t v : set) {
        profit_sum += prob.profit(v);
    }

    double path_cost = held_karp(source, target, set, all_pairs, n);
    return path_cost + profit_sum;
}

/// Try to grow an infeasible set by adding high-y candidates.
/// Returns the grown set (still infeasible, lb > ub).
void greedy_grow(const Problem& prob, std::span<const double> all_pairs, double ub,
                 std::vector<int32_t>& set, const std::vector<int32_t>& candidates,
                 int32_t max_dp_size = 15) {
    for (int v : candidates) {
        if (static_cast<int32_t>(set.size()) >= max_dp_size) {
            break;
        }

        // Skip if already in set
        bool in_set = false;
        for (int32_t u : set) {
            if (u == v) {
                in_set = true;
                break;
            }
        }
        if (in_set) {
            continue;
        }

        // Try adding
        set.push_back(v);
        double lb = compute_set_lb(prob, all_pairs, set, max_dp_size);
        if (lb <= ub + 1e-6) {
            set.pop_back();  // didn't help — remove
        }
        // else: kept, set is still infeasible
    }
}

/// Shrink an infeasible set to a minimal infeasible subset.
/// Repeatedly removes the node whose removal keeps lb highest above ub,
/// preferring to remove low-y nodes (which increases cut violation).
/// Will not shrink below min_size.
void greedy_shrink(const Problem& prob, std::span<const double> all_pairs, double ub,
                   std::vector<int32_t>& set, std::span<const double> y_values,
                   size_t min_size = 2) {
    bool improved = true;
    while (improved && set.size() > min_size) {
        improved = false;
        double best_lb = -kInf;
        size_t best_idx = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            // Try removing set[i]
            std::vector<int32_t> reduced;
            reduced.reserve(set.size() - 1);
            for (size_t j = 0; j < set.size(); ++j) {
                if (j != i) {
                    reduced.push_back(set[j]);
                }
            }

            double lb = compute_set_lb(prob, all_pairs, reduced);
            if (lb > ub + 1e-6) {
                // Prefer removing the node that:
                //  1) keeps lb highest (most slack for further shrinking)
                //  2) among ties, has the lowest y (maximizes violation)
                if (lb > best_lb + 1e-9 ||
                    (lb > best_lb - 1e-9 && y_values[set[i]] < y_values[set[best_idx]])) {
                    best_lb = lb;
                    best_idx = i;
                    improved = true;
                }
            }
        }

        if (improved) {
            set.erase(set.begin() + static_cast<ptrdiff_t>(best_idx));
        }
    }
}

/// Lift an infeasible pair S = {a, b} into a clique of pairwise conflicts.
///
/// The cut Σ_{v∈C} y_v ≤ 1 over the support C is valid exactly when every
/// two members of C form an infeasible set: a route visiting two of them
/// visits such a pair, and lb of a subset bounds the route from below.
/// Lifting is therefore sequential: j joins C only if it conflicts with every
/// node already in C, previously lifted nodes included (#171).  Candidates
/// come in `order` (y descending), so a high-y node is not blocked by a
/// low-y one lifted before it.
///
/// Covers of three or more nodes are emitted unlifted: the same rule there
/// checks every (|S|-1)-subset of C, whose count grows exponentially.
///
/// Returns lifted indices and values (coefficients).  RHS stays 1.
void lift_pair_cut(const Problem& prob, std::span<const double> all_pairs, double ub,
                   const std::vector<int32_t>& pair, const std::vector<int32_t>& order,
                   std::vector<int32_t>& out_nodes, std::vector<double>& out_coeffs) {
    out_nodes = pair;
    out_coeffs.assign(pair.size(), 1.0);

    for (int32_t j : order) {
        if (std::ranges::find(pair, j) != pair.end()) {
            continue;
        }
        const bool conflicts_with_all = std::ranges::all_of(out_nodes, [&](int32_t v) {
            return compute_set_lb(prob, all_pairs, {v, j}) > ub + 1e-6;
        });
        if (conflicts_with_all) {
            out_nodes.push_back(j);
            out_coeffs.push_back(1.0);
        }
    }
}

}  // namespace

std::vector<Cut> SPISeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double ub = ctx.upper_bound;

    // Need both all-pairs matrix and a finite upper bound
    if (ub >= kInf || ctx.all_pairs.empty()) {
        return {};
    }
    if (static_cast<int32_t>(ctx.all_pairs.size()) < n * n) {
        return {};
    }

    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const bool is_tour = prob.is_tour();

    // ── Collect candidates sorted by y descending ──
    // Non-terminals with positive y.  Sorting by y focuses the search on
    // nodes the LP most wants to select — proving their mutual infeasibility
    // gives the most useful cuts.
    struct Candidate {
        int32_t node;
        double y;
    };
    std::vector<Candidate> cands;
    for (int32_t i = 0; i < n; ++i) {
        if (i == source) {
            continue;
        }
        if (!is_tour && i == target) {
            continue;
        }
        if (ctx.y_values[i] <= tol) {
            continue;
        }
        cands.push_back({i, ctx.y_values[i]});
    }
    std::ranges::sort(cands, [](const Candidate& a, const Candidate& b) { return a.y > b.y; });

    // Flat candidate node list (y-descending) for grow phase
    const auto grow_cap = static_cast<size_t>(max_grow_candidates_);
    std::vector<int32_t> cand_nodes;
    cand_nodes.reserve(std::min(cands.size(), grow_cap));
    for (size_t i = 0; i < std::min(cands.size(), grow_cap); ++i) {
        cand_nodes.push_back(cands[i].node);
    }

    // Every candidate node (y-descending) for pair lifting
    std::vector<int32_t> lift_order;
    lift_order.reserve(cands.size());
    for (const auto& c : cands) {
        lift_order.push_back(c.node);
    }

    // ── Phase 1: Find all violated infeasible pairs ──
    std::vector<Cut> cuts;
    std::vector<std::vector<int32_t>> emitted;  // sorted supports, both phases
    auto is_duplicate = [&](const std::vector<int32_t>& set) -> bool {
        return std::ranges::any_of(emitted, [&](const auto& s) { return s == set; });
    };

    // Track infeasible pairs for seeding Phase 2
    struct InfeasiblePair {
        int32_t i, j;
        double violation;
    };
    std::vector<InfeasiblePair> seeds;

    for (size_t ai = 0; ai < cands.size(); ++ai) {
        const int32_t i = cands[ai].node;
        const double yi = cands[ai].y;

        for (size_t aj = ai + 1; aj < cands.size(); ++aj) {
            const int32_t j = cands[aj].node;
            const double yj = cands[aj].y;

            double violation = yi + yj - 1.0;
            if (violation <= tol) {
                continue;
            }

            double lb = compute_set_lb(prob, ctx.all_pairs, {i, j}, max_dp_size_);
            if (lb > ub + 1e-6) {
                // Lift: add nodes that conflict with every node in the cut
                std::vector<int32_t> lifted_nodes;
                std::vector<double> lifted_coeffs;
                lift_pair_cut(prob, ctx.all_pairs, ub, {i, j}, lift_order, lifted_nodes,
                              lifted_coeffs);
                seeds.push_back({i, j, violation});

                // Every violated pair of one clique lifts into the same row
                // (#176); emit it once.  Lifting still runs per pair: the
                // greedy clique grown from {i, j} need not be one emitted.
                auto sorted_support = lifted_nodes;
                std::ranges::sort(sorted_support);
                if (is_duplicate(sorted_support)) {
                    continue;
                }

                double sum_y = 0.0;
                Cut cut;
                cut.rhs = 1.0;  // |S|-1 where S is the minimal set (pair)
                for (size_t li = 0; li < lifted_nodes.size(); ++li) {
                    cut.indices.push_back(ctx.y_offset + lifted_nodes[li]);
                    cut.values.push_back(lifted_coeffs[li]);
                    sum_y += lifted_coeffs[li] * ctx.y_values[lifted_nodes[li]];
                }
                cut.violation = sum_y - cut.rhs;
                if (cut.violation > tol) {
                    emitted.push_back(std::move(sorted_support));
                    cuts.push_back(std::move(cut));
                }
            }
        }
    }

    // ── Phase 2: Greedy grow + shrink from infeasible pairs ──
    //
    // For each infeasible pair, greedily extend with high-y candidates,
    // then shrink to a minimal infeasible subset.
    //
    // Why this works: violation = Σy_v - (|S|-1).  For fractional nodes,
    // removing any v increases violation by (1-y_v) > 0.  So the minimal
    // infeasible set gives the strongest (most violated) cut.

    // Sort seeds by violation descending — most violated pairs first
    std::ranges::sort(seeds, [](const InfeasiblePair& a, const InfeasiblePair& b) {
        return a.violation > b.violation;
    });

    size_t num_seeds = std::min(seeds.size(), static_cast<size_t>(max_seeds_));
    for (size_t si = 0; si < num_seeds; ++si) {
        auto [pi, pj, _] = seeds[si];

        std::vector<int32_t> set = {pi, pj};

        // Grow: add high-y candidates that maintain infeasibility
        greedy_grow(prob, ctx.all_pairs, ub, set, cand_nodes, max_dp_size_);

        if (static_cast<int32_t>(set.size()) <= 2) {
            continue;  // pair covered by its lifted Phase 1 clique row
        }

        // Shrink: remove nodes to find minimal infeasible subset.
        // min_size=3: don't shrink back to pairs (each is covered by its
        // lifted Phase 1 clique row).
        greedy_shrink(prob, ctx.all_pairs, ub, set, ctx.y_values, /*min_size=*/3);

        // Compute violation
        double sum_y = 0.0;
        for (int32_t v : set) {
            sum_y += ctx.y_values[v];
        }
        double violation = sum_y - (static_cast<double>(set.size()) - 1.0);
        if (violation <= tol) {
            continue;
        }

        // Dedup
        auto sorted_set = set;
        std::ranges::sort(sorted_set);
        if (is_duplicate(sorted_set)) {
            continue;
        }
        emitted.push_back(sorted_set);

        // Emit the unlifted cut Σ_{v∈S} y_v ≤ |S| - 1 (see lift_pair_cut)
        Cut cut;
        cut.rhs = static_cast<double>(set.size()) - 1.0;
        for (int32_t v : set) {
            cut.indices.push_back(ctx.y_offset + v);
            cut.values.push_back(1.0);
        }
        cut.violation = violation;
        cuts.push_back(std::move(cut));
    }

    return cuts;
}

}  // namespace cptp::sep
