#include "sep/rglm_separator.h"

#include "core/gomory_hu.h"
#include "core/problem.h"
#include "sep/capacity_rows.h"

#include <cmath>
#include <span>
#include <vector>

namespace cptp::sep {

std::vector<Cut> RGLMSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.source();
    const bool is_tour = prob.is_tour();
    const int32_t path_target = prob.target();

    // Rows that can give a cut: 0 < Q_k < kInfiniteCapacity, as row 0 alone
    // was before, plus the surrogate row when it is on (#112).
    const SurrogateRow surrogate = surrogate_ ? surrogate_row(prob) : SurrogateRow{};
    const std::vector<CapacityRowView> rows = capacity_row_views(prob, rows_, surrogate);
    if (rows.empty()) {
        return {};
    }

    // Precompute each row's total demand (constant across all targets).
    std::vector<double> d_totals(rows.size(), 0.0);
    for (size_t r = 0; r < rows.size(); ++r) {
        for (int32_t i = 0; i < n; ++i) {
            if (i != depot && (is_tour || i != path_target)) {
                d_totals[r] += rows[r].demand[i];
            }
        }
    }

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || (!is_tour && target == path_target) || ctx.y_values[target] <= tol) {
            continue;
        }

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);

        // One cut per row (kAll), or the most violated row's (otherwise).
        Cut best;
        bool have_best = false;

        for (size_t row_idx = 0; row_idx < rows.size(); ++row_idx) {
            const std::span<const double> d = rows[row_idx].demand;
            const double Q = rows[row_idx].capacity;
            const double d_total = d_totals[row_idx];

            // S = {u : !reachable[u]} (target side, not containing depot)
            // Compute d_S = sum of row demands in S.
            double d_S = 0.0;
            for (int32_t i = 0; i < n; ++i) {
                if (!reachable[i] && (is_tour || i != path_target)) {
                    d_S += d[i];
                }
            }

            // alpha = d_S + 2*(d_total - d_S) = 2*d_total - d_S
            double alpha = 2.0 * d_total - d_S;
            double r = std::fmod(alpha, Q);
            if (r <= tol) {
                continue;  // no rounding benefit
            }

            double k = std::ceil(alpha / Q);
            if (k <= 1.0) {
                continue;  // no strengthening
            }

            // Compute cut_flow and star_sum over crossing edges.
            double cut_flow = 0.0;
            double star_sum = 0.0;  // sum of d_out(e)*x_e for e in delta(S)

            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                bool u_in_S = !reachable[u] && (is_tour || u != path_target);
                bool v_in_S = !reachable[v] && (is_tour || v != path_target);
                if (u_in_S == v_in_S) {
                    continue;
                }

                int32_t outside = u_in_S ? v : u;
                double d_out = d[outside];

                cut_flow += ctx.x_values[e];
                star_sum += d_out * ctx.x_values[e];
            }

            // beta = sum_{i in S} d_i*(1 - y_i) + 2*(d_total - d_S) - star_sum
            double d_NmS = d_total - d_S;
            double beta_y_part = 0.0;
            for (int32_t i = 0; i < n; ++i) {
                if (!reachable[i] && (is_tour || i != path_target)) {
                    beta_y_part += d[i] * (1.0 - ctx.y_values[i]);
                }
            }
            double beta = beta_y_part + 2.0 * d_NmS - star_sum;

            // Violation: (2k - 2*beta/r) - cut_flow > tol
            double violation = (2.0 * k - 2.0 * beta / r) - cut_flow;
            if (violation <= tol) {
                continue;
            }

            // Build cut in <= form:
            //   (2/r)*sum_{i in S} d_i*y_i
            //     + sum_{e in delta(S)} (2*d_out/r - 1)*x_e
            //     <= 2*(alpha/r - k)
            Cut cut;

            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                bool u_in_S = !reachable[u] && (is_tour || u != path_target);
                bool v_in_S = !reachable[v] && (is_tour || v != path_target);
                if (u_in_S == v_in_S) {
                    continue;
                }

                int32_t outside = u_in_S ? v : u;
                double d_out = d[outside];

                double coeff = 2.0 * d_out / r - 1.0;
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(coeff);
            }

            for (int32_t i = 0; i < n; ++i) {
                if (reachable[i] || (!is_tour && i == path_target)) {
                    continue;
                }
                double d_i = d[i];
                if (d_i <= 0) {
                    continue;
                }

                double coeff = 2.0 * d_i / r;
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(coeff);
            }

            cut.rhs = 2.0 * (alpha / r - k);
            cut.violation = violation;
            if (cut.indices.empty()) {
                continue;
            }
            if (rows_ != CapacityRows::kMostViolated) {
                all_cuts.push_back(std::move(cut));
            } else if (!have_best || cut.violation > best.violation) {
                best = std::move(cut);
                have_best = true;
            }
        }

        if (have_best) {
            all_cuts.push_back(std::move(best));
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
