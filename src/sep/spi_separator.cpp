#include "sep/spi_separator.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "core/problem.h"

namespace rcspp::sep {

std::vector<Cut> SPISeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double ub = ctx.upper_bound;
    constexpr double inf = std::numeric_limits<double>::infinity();

    // Need both all-pairs matrix and a finite upper bound
    if (ub >= inf || ctx.all_pairs.empty()) return {};
    if (static_cast<int32_t>(ctx.all_pairs.size()) < n * n) return {};

    const auto& d = ctx.all_pairs;
    const int32_t source = prob.source();
    const int32_t target = prob.target();
    const bool is_tour = prob.is_tour();
    const double correction = is_tour ? prob.profit(source) : 0.0;

    std::vector<Cut> cuts;

    // Collect candidate nodes: non-terminals with positive y value.
    // Terminals are always visited, so pair (terminal, i) reduces to
    // node elimination (already handled by edge elimination preprocessing).
    std::vector<int32_t> candidates;
    for (int32_t i = 0; i < n; ++i) {
        if (i == source) continue;
        if (!is_tour && i == target) continue;
        if (ctx.y_values[i] <= tol) continue;
        candidates.push_back(i);
    }

    // ── Pair inequalities: y_i + y_j <= 1 ──
    //
    // For a tour visiting depot (d), i, and j, the route decomposes into
    // three segments: d→i, i→j, j→d.  Each segment's shortest-path cost
    // d[u][v] includes the profits of all nodes u,...,v.  Junction nodes
    // (d, i, j) appear in two segments, so their profit is subtracted twice.
    // Adding it back once gives the correct tour cost:
    //   actual_obj = d[d][i] + d[i][j] + d[j][d] + profit(d) + profit(i) + profit(j)
    //
    // For a path (source!=target), terminals appear only once (at the ends),
    // so only the interior junction nodes (i, j) need the correction:
    //   actual_obj = d[s][i] + d[i][j] + d[j][t] + profit(i) + profit(j)

    for (size_t ai = 0; ai < candidates.size(); ++ai) {
        const int32_t i = candidates[ai];
        const double yi = ctx.y_values[i];
        const double pi = prob.profit(i);

        for (size_t aj = ai + 1; aj < candidates.size(); ++aj) {
            const int32_t j = candidates[aj];
            const double yj = ctx.y_values[j];

            double violation = yi + yj - 1.0;
            if (violation <= tol) continue;

            const double pj = prob.profit(j);

            double lb;
            if (is_tour) {
                double base = correction + pi + pj;
                double lb1 = d[source * n + i] + d[i * n + j] + d[j * n + source] + base;
                double lb2 = d[source * n + j] + d[j * n + i] + d[i * n + source] + base;
                lb = std::min(lb1, lb2);
            } else {
                double base = pi + pj;
                double lb1 = d[source * n + i] + d[i * n + j] + d[j * n + target] + base;
                double lb2 = d[source * n + j] + d[j * n + i] + d[i * n + target] + base;
                lb = std::min(lb1, lb2);
            }

            if (lb > ub + 1e-6) {
                Cut cut;
                cut.violation = violation;
                cut.indices = {ctx.y_offset + i, ctx.y_offset + j};
                cut.values = {1.0, 1.0};
                cut.rhs = 1.0;
                cuts.push_back(std::move(cut));
            }
        }
    }

    // ── Triplet inequalities: y_i + y_j + y_k <= 2 ──
    //
    // Same decomposition with 4 segments and 3 interior junction nodes.
    // Tour: lb = d[d][i] + d[i][j] + d[j][k] + d[k][d] + correction + pi + pj + pk
    // Path: lb = d[s][i] + d[i][j] + d[j][k] + d[k][t] + pi + pj + pk
    //
    // We enumerate all 3 cyclic orderings (tour) or 6 permutations (path)
    // and take the minimum.  Skipped for large candidate sets to avoid O(n³).
    constexpr size_t kMaxTripletCandidates = 200;

    if (candidates.size() <= kMaxTripletCandidates) {
        for (size_t ai = 0; ai < candidates.size(); ++ai) {
            const int32_t i = candidates[ai];
            const double yi = ctx.y_values[i];
            const double pi = prob.profit(i);

            for (size_t aj = ai + 1; aj < candidates.size(); ++aj) {
                const int32_t j = candidates[aj];
                const double yj = ctx.y_values[j];

                // Early exit: if y_i + y_j + 1.0 <= 2 + tol, no triple can violate
                if (yi + yj <= 1.0 + tol) continue;

                const double pj = prob.profit(j);

                for (size_t ak = aj + 1; ak < candidates.size(); ++ak) {
                    const int32_t k = candidates[ak];
                    const double yk = ctx.y_values[k];

                    double violation = yi + yj + yk - 2.0;
                    if (violation <= tol) continue;

                    const double pk = prob.profit(k);

                    double lb;
                    if (is_tour) {
                        double base = correction + pi + pj + pk;
                        // 3 distinct cyclic orderings (i-j-k, i-k-j, j-i-k)
                        double lb1 = d[source*n+i] + d[i*n+j] + d[j*n+k] + d[k*n+source] + base;
                        double lb2 = d[source*n+i] + d[i*n+k] + d[k*n+j] + d[j*n+source] + base;
                        double lb3 = d[source*n+j] + d[j*n+i] + d[i*n+k] + d[k*n+source] + base;
                        lb = std::min({lb1, lb2, lb3});
                    } else {
                        double base = pi + pj + pk;
                        // 6 permutations of {i,j,k}
                        double lb1 = d[source*n+i] + d[i*n+j] + d[j*n+k] + d[k*n+target] + base;
                        double lb2 = d[source*n+i] + d[i*n+k] + d[k*n+j] + d[j*n+target] + base;
                        double lb3 = d[source*n+j] + d[j*n+i] + d[i*n+k] + d[k*n+target] + base;
                        double lb4 = d[source*n+j] + d[j*n+k] + d[k*n+i] + d[i*n+target] + base;
                        double lb5 = d[source*n+k] + d[k*n+i] + d[i*n+j] + d[j*n+target] + base;
                        double lb6 = d[source*n+k] + d[k*n+j] + d[j*n+i] + d[i*n+target] + base;
                        lb = std::min({lb1, lb2, lb3, lb4, lb5, lb6});
                    }

                    if (lb > ub + 1e-6) {
                        Cut cut;
                        cut.violation = violation;
                        cut.indices = {ctx.y_offset + i, ctx.y_offset + j, ctx.y_offset + k};
                        cut.values = {1.0, 1.0, 1.0};
                        cut.rhs = 2.0;
                        cuts.push_back(std::move(cut));
                    }
                }
            }
        }
    }

    return cuts;
}

}  // namespace rcspp::sep
