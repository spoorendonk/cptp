#include "sep/multistar_separator.h"

#include "core/gomory_hu.h"
#include "core/problem.h"

#include <cmath>
#include <vector>

namespace cptp::sep {

std::vector<Cut> MultistarSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.source();
    const bool is_tour = prob.is_tour();
    const int32_t path_target = prob.target();

    // Rows that can give a cut: 0 < Q_k < kInfiniteCapacity, as row 0 alone
    // was before.
    std::vector<int32_t> rows;
    const int32_t num_rows = rows_ == CapacityRows::kFirst ? 1 : prob.num_capacities();
    for (int32_t k = 0; k < num_rows; ++k) {
        const double Q = prob.capacity(k);
        if (Q > 0 && Q < kInfiniteCapacity) {
            rows.push_back(k);
        }
    }
    if (rows.empty()) {
        return {};
    }

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || (!is_tour && target == path_target) || ctx.y_values[target] <= tol) {
            continue;
        }

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);

        // S = {u : !reachable[u]} (target side, not containing depot)
        // GLM for CPTP, per row k (d = d^k, Q = Q_k):
        //   sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e
        //     - (2/Q)*sum_{i in S} q_i*y_i >= 0
        // In <= form:
        //   (2/Q)*sum_{i in S} q_i*y_i
        //     - sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e <= 0
        // One cut per row (kAll), or the most violated row's (otherwise).
        Cut best;
        bool have_best = false;

        for (const int32_t row : rows) {
            const double Q = prob.capacity(row);
            double lhs = 0.0;
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
                double q_outside = prob.demand(row, outside);

                double coeff = 1.0 - 2.0 * q_outside / Q;
                lhs -= coeff * ctx.x_values[e];
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-coeff);
            }

            for (int32_t i = 0; i < n; ++i) {
                if (reachable[i] || (!is_tour && i == path_target)) {
                    continue;
                }
                double q_i = prob.demand(row, i);
                if (q_i <= 0) {
                    continue;
                }

                double coeff = 2.0 * q_i / Q;
                lhs += coeff * ctx.y_values[i];
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(coeff);
            }

            cut.rhs = 0.0;

            if (lhs <= tol || cut.indices.empty()) {
                continue;
            }
            cut.violation = lhs;
            if (rows_ == CapacityRows::kAll) {
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
