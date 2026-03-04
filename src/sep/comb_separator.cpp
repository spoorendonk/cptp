#include "sep/comb_separator.h"

#include <algorithm>
#include <vector>

#include "core/problem.h"

namespace cptp::sep {

std::vector<Cut> CombSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.source();
    const bool is_tour = prob.is_tour();
    const int32_t path_target = prob.target();

    std::vector<Cut> cuts;

    // Try two BFS thresholds for handle construction.
    for (double threshold : {0.5, 0.3}) {
        // BFS from depot to find handle set H.
        std::vector<bool> in_handle(n, false);
        std::vector<int32_t> queue;
        queue.push_back(depot);
        in_handle[depot] = true;
        if (!is_tour && path_target != depot) {
            // In s-t mode, treat the split depot terminals as one root.
            queue.push_back(path_target);
            in_handle[path_target] = true;
        }

        for (size_t qi = 0; qi < queue.size(); ++qi) {
            int32_t u = queue[qi];
            for (auto e : graph.incident_edges(u)) {
                int32_t v = graph.other_endpoint(e, u);
                if (!in_handle[v] && ctx.x_values[e] > threshold) {
                    in_handle[v] = true;
                    queue.push_back(v);
                }
            }
        }

        int32_t handle_size = static_cast<int32_t>(queue.size());
        if (handle_size < 2 || handle_size >= n) continue;

        // Find candidate teeth: edges crossing the handle boundary.
        // Each tooth T_i = {u_i, v_i} with u_i in H, v_i not in H.
        // Contribution to violation: x_e - y_{v_i}.
        struct Tooth {
            int32_t inside;
            int32_t outside;
            int32_t edge;
            double contrib;
        };
        std::vector<Tooth> teeth;

        for (auto e : graph.edges()) {
            double xval = ctx.x_values[e];
            if (xval < tol) continue;

            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);

            int32_t ins = -1, out = -1;
            if (in_handle[u] && !in_handle[v]) {
                ins = u;
                out = v;
            } else if (in_handle[v] && !in_handle[u]) {
                ins = v;
                out = u;
            } else {
                continue;
            }
            if (!is_tour && out == path_target) continue;

            teeth.push_back({ins, out, e, xval - ctx.y_values[out]});
        }

        // Sort by contribution descending (most violating first).
        // stable_sort preserves deterministic edge-iteration order for equal contributions.
        std::stable_sort(teeth.begin(), teeth.end(),
                         [](const Tooth& a, const Tooth& b) {
                             return a.contrib > b.contrib;
                         });

        // Greedy selection: no shared inside or outside nodes.
        std::vector<bool> used_inside(n, false);
        std::vector<bool> used_outside(n, false);
        std::vector<Tooth> selected;

        auto inside_key = [&](int32_t node) {
            if (!is_tour && node == path_target) return depot;
            return node;
        };

        for (const auto& tooth : teeth) {
            const int32_t in_key = inside_key(tooth.inside);
            if (used_inside[in_key] || used_outside[tooth.outside])
                continue;
            selected.push_back(tooth);
            used_inside[in_key] = true;
            used_outside[tooth.outside] = true;
        }

        // Need >= 3 teeth, odd count.
        if (selected.size() < 3) continue;
        size_t t = selected.size();
        if (t % 2 == 0) t--;
        selected.resize(t);

        // Evaluate: x(E(H)) + Σ x_{tooth_i} - Σ_{j∈H} y_j - Σ y_{outside_i}
        double lhs = 0.0;

        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            if (!is_tour &&
                ((u == depot && v == path_target) ||
                 (u == path_target && v == depot))) {
                continue;
            }
            if (in_handle[u] && in_handle[v]) {
                lhs += ctx.x_values[e];
            }
        }

        for (const auto& tooth : selected) {
            lhs += ctx.x_values[tooth.edge];
        }

        for (int32_t i = 0; i < n; ++i) {
            if (!is_tour && i == path_target) continue;
            if (in_handle[i]) lhs -= ctx.y_values[i];
        }

        for (const auto& tooth : selected) {
            lhs -= ctx.y_values[tooth.outside];
        }

        double rhs = static_cast<double>(t - 1) / 2.0;

        if (lhs - rhs <= tol) continue;

        // Emit cut: a^T x <= rhs.
        Cut cut;
        cut.violation = lhs - rhs;

        // E(H) edges: coeff +1.
        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            if (!is_tour &&
                ((u == depot && v == path_target) ||
                 (u == path_target && v == depot))) {
                continue;
            }
            if (in_handle[u] && in_handle[v]) {
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(1.0);
            }
        }

        // Tooth edges: coeff +1.
        for (const auto& tooth : selected) {
            cut.indices.push_back(ctx.x_offset + tooth.edge);
            cut.values.push_back(1.0);
        }

        // y_j for j in H: coeff -1.
        for (int32_t i = 0; i < n; ++i) {
            if (!is_tour && i == path_target) continue;
            if (in_handle[i]) {
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(-1.0);
            }
        }

        // y_{outside_i} for each tooth: coeff -1.
        for (const auto& tooth : selected) {
            cut.indices.push_back(ctx.y_offset + tooth.outside);
            cut.values.push_back(-1.0);
        }

        cut.rhs = rhs;
        cuts.push_back(std::move(cut));
    }

    return cuts;
}

}  // namespace cptp::sep
