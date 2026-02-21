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
        double sec_violation = 2.0 * yi - flow;
        if (sec_violation <= tol) continue;

        // Count |E(S)| and |δ(S)| to choose sparser form.
        // Also count |S| for the inside form's y-variable terms.
        int32_t n_inside = 0;
        int32_t n_delta = 0;
        int32_t s_size = 0;
        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            bool u_in = !reachable[u];
            bool v_in = !reachable[v];
            if (u_in && v_in) n_inside++;
            else if (u_in != v_in) n_delta++;
        }
        for (int32_t i = 0; i < n; ++i) {
            if (!reachable[i]) s_size++;
        }

        // Cut form:   |δ(S)| + 1 nonzeros
        // Inside form: |E(S)| + |S| - 1 nonzeros
        bool use_inside = (n_inside + s_size - 1) < (n_delta + 1);

        Cut cut;
        cut.violation = sec_violation;

        if (use_inside) {
            // Inside form (via degree substitution x(δ({i}))=2y_i):
            //   x(E(S)) - Σ_{j∈S\{target}} y_j  ≤  0
            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                if (!reachable[u] && !reachable[v]) {
                    cut.indices.push_back(ctx.x_offset + e);
                    cut.values.push_back(1.0);
                }
            }
            for (int32_t i = 0; i < n; ++i) {
                if (reachable[i] || i == target) continue;
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(-1.0);
            }
            cut.rhs = 0.0;
        } else {
            // Cut form: 2·y_target - x(δ(S)) ≤ 0
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
        }

        if (!cut.indices.empty()) {
            all_cuts.push_back(std::move(cut));
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
