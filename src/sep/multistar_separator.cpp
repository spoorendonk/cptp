#include "sep/multistar_separator.h"

#include <cmath>

#include "core/gomory_hu.h"
#include "core/problem.h"

namespace cptp::sep {

std::vector<Cut> MultistarSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();

    if (Q <= 0 || Q >= 1e17) return {};

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || ctx.y_values[target] <= tol) continue;

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);

        // S = {u : !reachable[u]} (target side, not containing depot)
        // GLM for CPTP:
        //   sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e
        //     - (2/Q)*sum_{i in S} q_i*y_i >= 0
        // In <= form:
        //   (2/Q)*sum_{i in S} q_i*y_i
        //     - sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e <= 0

        double lhs = 0.0;
        Cut cut;

        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            bool u_in_S = !reachable[u];
            bool v_in_S = !reachable[v];
            if (u_in_S == v_in_S) continue;

            int32_t outside = u_in_S ? v : u;
            double q_outside = prob.demand(outside);

            double coeff = 1.0 - 2.0 * q_outside / Q;
            lhs -= coeff * ctx.x_values[e];
            cut.indices.push_back(ctx.x_offset + e);
            cut.values.push_back(-coeff);
        }

        for (int32_t i = 0; i < n; ++i) {
            if (reachable[i]) continue;
            double q_i = prob.demand(i);
            if (q_i <= 0) continue;

            double coeff = 2.0 * q_i / Q;
            lhs += coeff * ctx.y_values[i];
            cut.indices.push_back(ctx.y_offset + i);
            cut.values.push_back(coeff);
        }

        cut.rhs = 0.0;

        if (lhs > tol && !cut.indices.empty()) {
            cut.violation = lhs;
            all_cuts.push_back(std::move(cut));
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
