#include "sep/spi_separator.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "core/problem.h"

namespace cptp::sep {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Default constants moved to SPISeparator members; these are used as
// fallback defaults in the anonymous-namespace helpers that take them as params.

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

double held_karp(int32_t source, int32_t target,
                 const std::vector<int32_t>& nodes,
                 std::span<const double> ap, int32_t n) {
    const int32_t k = static_cast<int32_t>(nodes.size());

    if (k == 0) return ap[static_cast<size_t>(source) * n + target];

    auto d = [&](int32_t u, int32_t v) -> double {
        return ap[static_cast<size_t>(u) * n + v];
    };

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
            if (!(mask & (1 << last))) continue;
            double cur = dp[static_cast<size_t>(mask) * k + last];
            if (cur >= kInf) continue;

            for (int32_t next = 0; next < k; ++next) {
                if (mask & (1 << next)) continue;
                int32_t new_mask = mask | (1 << next);
                double cost = cur + d(nodes[last], nodes[next]);
                auto& slot = dp[static_cast<size_t>(new_mask) * k + next];
                if (cost < slot) slot = cost;
            }
        }
    }

    // Complete: find best way to reach target
    const int32_t full = states - 1;
    double best = kInf;
    for (int32_t last = 0; last < k; ++last) {
        double cost = dp[static_cast<size_t>(full) * k + last] + d(nodes[last], target);
        if (cost < best) best = cost;
    }
    return best;
}

/// Compute the lower bound on the objective of any tour/path visiting all
/// nodes in `set`.  Returns -∞ if the set is too large for exact DP.
double compute_set_lb(const Problem& prob,
                      std::span<const double> all_pairs,
                      const std::vector<int32_t>& set,
                      int32_t max_dp_size = 15) {
    const int32_t n = prob.num_nodes();
    const int32_t k = static_cast<int32_t>(set.size());
    if (k == 0 || k > max_dp_size) return -kInf;

    const int32_t source = prob.source();
    const int32_t target = prob.is_tour() ? source : prob.target();
    const double correction = prob.is_tour() ? prob.profit(source) : 0.0;

    // Junction-node profit correction: each node in the set appears in two
    // consecutive segments, so its profit is subtracted twice.  Add it back.
    double profit_sum = correction;
    for (int32_t v : set) profit_sum += prob.profit(v);

    double path_cost = held_karp(source, target, set, all_pairs, n);
    return path_cost + profit_sum;
}

/// Try to grow an infeasible set by adding high-y candidates.
/// Returns the grown set (still infeasible, lb > ub).
void greedy_grow(const Problem& prob,
                 std::span<const double> all_pairs,
                 double ub,
                 std::vector<int32_t>& set,
                 const std::vector<int32_t>& candidates,
                 int32_t max_dp_size = 15) {
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        if (static_cast<int32_t>(set.size()) >= max_dp_size) break;

        int32_t v = candidates[ci];

        // Skip if already in set
        bool in_set = false;
        for (int32_t u : set) {
            if (u == v) { in_set = true; break; }
        }
        if (in_set) continue;

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
void greedy_shrink(const Problem& prob,
                   std::span<const double> all_pairs,
                   double ub,
                   std::vector<int32_t>& set,
                   std::span<const double> y_values,
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
                if (j != i) reduced.push_back(set[j]);
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

/// Lift a minimal infeasible set S by scanning nodes j ∉ S.
/// Node j gets coefficient 1 if replacing every i ∈ S with j still yields
/// an infeasible set (lb > ub).  This is the "extended cover" lifting:
/// any node that is at least as "expensive" as every member of the minimal
/// cover can be added with coefficient 1.
///
/// Returns lifted indices and values (coefficients).  RHS stays |S|-1.
void lift_cut(const Problem& prob,
              std::span<const double> all_pairs,
              double ub,
              const std::vector<int32_t>& set,
              std::span<const double> y_values,
              int32_t source, int32_t target, bool is_tour,
              double tol,
              std::vector<int32_t>& out_nodes,
              std::vector<double>& out_coeffs) {
    const int32_t n = prob.num_nodes();

    // Start with the core set (coefficient 1)
    out_nodes = set;
    out_coeffs.assign(set.size(), 1.0);

    // For each candidate j ∉ S, check all |S| replacements
    for (int32_t j = 0; j < n; ++j) {
        if (j == source) continue;
        if (!is_tour && j == target) continue;
        if (y_values[j] <= tol) continue;

        // Skip if j ∈ S
        bool in_set = false;
        for (int32_t v : set) {
            if (v == j) { in_set = true; break; }
        }
        if (in_set) continue;

        // Check: for every i ∈ S, is (S \ {i}) ∪ {j} infeasible?
        bool can_lift = true;
        for (size_t idx = 0; idx < set.size(); ++idx) {
            std::vector<int32_t> swapped;
            swapped.reserve(set.size());
            for (size_t k = 0; k < set.size(); ++k) {
                swapped.push_back(k == idx ? j : set[k]);
            }
            double lb = compute_set_lb(prob, all_pairs, swapped);
            if (lb <= ub + 1e-6) {
                can_lift = false;
                break;
            }
        }

        if (can_lift) {
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
    if (ub >= kInf || ctx.all_pairs.empty()) return {};
    if (static_cast<int32_t>(ctx.all_pairs.size()) < n * n) return {};

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
        if (i == source) continue;
        if (!is_tour && i == target) continue;
        if (ctx.y_values[i] <= tol) continue;
        cands.push_back({i, ctx.y_values[i]});
    }
    std::sort(cands.begin(), cands.end(),
              [](const Candidate& a, const Candidate& b) { return a.y > b.y; });

    // Flat candidate node list (y-descending) for grow phase
    const size_t grow_cap = static_cast<size_t>(max_grow_candidates_);
    std::vector<int32_t> cand_nodes;
    cand_nodes.reserve(std::min(cands.size(), grow_cap));
    for (size_t i = 0; i < std::min(cands.size(), grow_cap); ++i) {
        cand_nodes.push_back(cands[i].node);
    }

    // ── Phase 1: Find all violated infeasible pairs ──
    std::vector<Cut> cuts;
    std::vector<std::vector<int32_t>> emitted;  // dedup

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
            if (violation <= tol) continue;

            double lb = compute_set_lb(prob, ctx.all_pairs, {i, j}, max_dp_size_);
            if (lb > ub + 1e-6) {
                // Lift: scan for nodes that can replace either member
                std::vector<int32_t> base_set = {i, j};
                std::vector<int32_t> lifted_nodes;
                std::vector<double> lifted_coeffs;
                lift_cut(prob, ctx.all_pairs, ub, base_set,
                         ctx.y_values, source, target, is_tour, tol,
                         lifted_nodes, lifted_coeffs);

                double sum_y = 0.0;
                Cut cut;
                cut.rhs = 1.0;  // |S|-1 where S is the minimal set (pair)
                for (size_t li = 0; li < lifted_nodes.size(); ++li) {
                    cut.indices.push_back(ctx.y_offset + lifted_nodes[li]);
                    cut.values.push_back(lifted_coeffs[li]);
                    sum_y += lifted_coeffs[li] * ctx.y_values[lifted_nodes[li]];
                }
                cut.violation = sum_y - cut.rhs;
                if (cut.violation > tol)
                    cuts.push_back(std::move(cut));

                seeds.push_back({i, j, violation});
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
    std::sort(seeds.begin(), seeds.end(),
              [](const InfeasiblePair& a, const InfeasiblePair& b) {
                  return a.violation > b.violation;
              });

    auto is_duplicate = [&](const std::vector<int32_t>& set) -> bool {
        for (const auto& s : emitted) {
            if (s == set) return true;
        }
        return false;
    };

    size_t num_seeds = std::min(seeds.size(), static_cast<size_t>(max_seeds_));
    for (size_t si = 0; si < num_seeds; ++si) {
        auto [pi, pj, _] = seeds[si];

        std::vector<int32_t> set = {pi, pj};

        // Grow: add high-y candidates that maintain infeasibility
        greedy_grow(prob, ctx.all_pairs, ub, set, cand_nodes, max_dp_size_);

        if (static_cast<int32_t>(set.size()) <= 2) continue;  // pair already emitted

        // Shrink: remove nodes to find minimal infeasible subset.
        // min_size=3: don't shrink back to pairs (already emitted in Phase 1).
        greedy_shrink(prob, ctx.all_pairs, ub, set, ctx.y_values, /*min_size=*/3);

        // Compute violation
        double sum_y = 0.0;
        for (int32_t v : set) sum_y += ctx.y_values[v];
        double violation = sum_y - (static_cast<double>(set.size()) - 1.0);
        if (violation <= tol) continue;

        // Dedup
        auto sorted_set = set;
        std::sort(sorted_set.begin(), sorted_set.end());
        if (is_duplicate(sorted_set)) continue;
        emitted.push_back(sorted_set);

        // Lift: scan for nodes that can replace any member of the minimal set
        std::vector<int32_t> lifted_nodes;
        std::vector<double> lifted_coeffs;
        lift_cut(prob, ctx.all_pairs, ub, set,
                 ctx.y_values, source, target, is_tour, tol,
                 lifted_nodes, lifted_coeffs);

        // Emit cut: Σ α_j y_j ≤ |S| - 1
        double lifted_sum_y = 0.0;
        Cut cut;
        cut.rhs = static_cast<double>(set.size()) - 1.0;
        for (size_t li = 0; li < lifted_nodes.size(); ++li) {
            cut.indices.push_back(ctx.y_offset + lifted_nodes[li]);
            cut.values.push_back(lifted_coeffs[li]);
            lifted_sum_y += lifted_coeffs[li] * ctx.y_values[lifted_nodes[li]];
        }
        cut.violation = lifted_sum_y - cut.rhs;
        if (cut.violation <= tol) continue;
        cuts.push_back(std::move(cut));
    }

    return cuts;
}

}  // namespace cptp::sep
