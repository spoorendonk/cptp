#include "sep/rci_separator.h"

#include <cmath>

#include "core/gomory_hu.h"
#include "core/problem.h"

namespace cptp::sep {

std::vector<Cut> RCISeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();

    if (Q <= 0) return {};

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || ctx.y_values[target] <= tol) continue;

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);

        // S = nodes NOT reachable from depot (target side)
        // Compute total demand of S
        double d_S = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            if (!reachable[i]) {
                d_S += prob.demand(i);
            }
        }

        double Q_r = std::fmod(d_S, Q);
        if (Q_r <= tol) continue;  // d(S) is multiple of Q

        double k = std::ceil(d_S / Q);
        if (k <= 1.0) continue;  // same as SEC

        // CPTP RCI: x(delta(S)) - (2/Q_r) * sum_{i in S} q_i * y_i
        //   >= 2*(ceil(d(S)/Q) - d(S)/Q_r)
        double served_demand = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            if (!reachable[i]) {
                served_demand += prob.demand(i) * ctx.y_values[i];
            }
        }

        double rhs = 2.0 * (k - d_S / Q_r);
        double lhs = flow - (2.0 * served_demand) / Q_r;

        double rci_violation = rhs - lhs;
        if (rci_violation <= tol) continue;  // not violated

        // In <= form: (2/Q_r)*sum q_i*y_i - x(delta(S)) <= -rhs
        Cut cut;
        cut.violation = rci_violation;

        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            if (reachable[u] != reachable[v]) {
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-1.0);
            }
        }

        for (int32_t i = 0; i < n; ++i) {
            if (reachable[i]) continue;
            double q_i = prob.demand(i);
            if (q_i <= 0) continue;

            double coeff = 2.0 * q_i / Q_r;
            cut.indices.push_back(ctx.y_offset + i);
            cut.values.push_back(coeff);
        }

        cut.rhs = -rhs;

        if (!cut.indices.empty()) {
            all_cuts.push_back(std::move(cut));
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
