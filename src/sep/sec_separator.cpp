#include "sep/sec_separator.h"

#include "core/gomory_hu.h"
#include "core/problem.h"

namespace rcspp::sep {

std::vector<Cut> SECSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t source = prob.source();
    const int32_t prob_target = prob.target();
    const bool is_tour = prob.is_tour();

    std::vector<Cut> all_cuts;

    for (int32_t target = 0; target < n; ++target) {
        if (target == source || ctx.y_values[target] <= tol) continue;

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);
        double yi = ctx.y_values[target];

        // S = non-root side (contains target). Root = source.
        // For paths: if S contains prob_target, the path can enter S and
        // terminate there, requiring only 1 cut crossing.
        // If S does NOT contain prob_target, path must enter and leave S → 2 crossings.
        bool target_in_S = !is_tour && !reachable[prob_target];
        double rhs_coeff = target_in_S ? 1.0 : 2.0;

        double sec_violation = rhs_coeff * yi - flow;
        if (sec_violation <= tol) continue;

        // Count |E(S)| and |δ(S)| to choose sparser form.
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
        // For paths with target_in_S, inside form derivation is complex,
        // so always use cut form in that case.
        // Inside form is degenerate (empty) when S = {target} (singleton),
        // so require at least 1 nonzero.
        int32_t inside_nnz = n_inside + s_size - 1;
        bool use_inside = !target_in_S && inside_nnz > 0 &&
                          inside_nnz < (n_delta + 1);

        Cut cut;
        cut.violation = sec_violation;

        if (use_inside) {
            // Inside form (via degree substitution x(δ({i}))=2y_i):
            //   x(E(S)) - Σ_{j∈S\{target}} y_j  ≤  0
            // Only used when S contains no path terminals (all deg=2).
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
            // Cut form: rhs_coeff·y_target - x(δ(S)) ≤ 0
            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                if (reachable[u] != reachable[v]) {
                    cut.indices.push_back(ctx.x_offset + e);
                    cut.values.push_back(-1.0);
                }
            }
            cut.indices.push_back(ctx.y_offset + target);
            cut.values.push_back(rhs_coeff);
            cut.rhs = 0.0;
        }

        if (!cut.indices.empty()) {
            all_cuts.push_back(std::move(cut));
        }
    }

    return all_cuts;
}

}  // namespace rcspp::sep
