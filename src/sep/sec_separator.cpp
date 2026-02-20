#include "sep/sec_separator.h"

#include "core/gomory_hu.h"
#include "core/problem.h"

namespace cptp::sep {

std::vector<Cut> SECSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.depot();

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || ctx.y_values[target] <= tol) continue;

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);
        double yi = ctx.y_values[target];

        // Undirected SEC: flow from depot to i >= 2*y_i
        if (flow >= 2.0 * yi - tol) continue;

        // SEC: sum(x_e for edges crossing cut) >= 2*y_target
        // In <= form: 2*y_target - sum(x_e) <= 0
        Cut cut;
        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            if (reachable[u] != reachable[v]) {
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-1.0);
            }
        }

        cut.indices.push_back(ctx.y_offset + target);
        cut.values.push_back(2.0);
        cut.rhs = 0.0;

        if (!cut.indices.empty()) {
            all_cuts.push_back(std::move(cut));
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
