#include "sep/comb_separator.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/problem.h"

namespace cptp::sep {

std::vector<Cut> CombSeparator::separate(const SeparationContext& ctx) {
    return {};  // DISABLED: comb cuts need validation for CPTP
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.depot();

    std::vector<Cut> cuts;

    // BFS from depot in the support graph to find the "handle" set
    std::vector<bool> in_handle(n, false);
    std::vector<int32_t> queue;
    queue.push_back(depot);
    in_handle[depot] = true;

    for (size_t qi = 0; qi < queue.size(); ++qi) {
        int32_t u = queue[qi];
        for (auto e : graph.incident_edges(u)) {
            int32_t v = graph.other_endpoint(e, u);
            if (!in_handle[v] && ctx.x_values[e] > 0.5 + tol) {
                in_handle[v] = true;
                queue.push_back(v);
            }
        }
    }

    // Find teeth: edges with fractional x crossing the handle boundary
    struct Tooth {
        int32_t inside;
        int32_t outside;
        int32_t edge;
        double x_val;
    };
    std::vector<Tooth> teeth;

    for (auto e : graph.edges()) {
        double xval = ctx.x_values[e];
        if (xval < tol || xval > 1.0 - tol) continue;

        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);

        if (in_handle[u] && !in_handle[v]) {
            teeth.push_back({u, v, e, xval});
        } else if (in_handle[v] && !in_handle[u]) {
            teeth.push_back({v, u, e, xval});
        }
    }

    // Sort teeth by violation potential (lowest x-value = most violated)
    std::sort(teeth.begin(), teeth.end(),
              [](const Tooth& a, const Tooth& b) { return a.x_val < b.x_val; });

    if (teeth.size() < 3) return cuts;

    // Select non-overlapping teeth
    std::vector<bool> used_outside(n, false);
    std::vector<Tooth> selected;

    for (const auto& tooth : teeth) {
        if (used_outside[tooth.outside]) continue;
        selected.push_back(tooth);
        used_outside[tooth.outside] = true;
    }

    if (selected.size() >= 3) {
        size_t k = selected.size();
        if (k % 2 == 0) k--;
        selected.resize(k);

        double lhs = 0.0;

        // x(delta(H)): edges crossing handle boundary
        for (auto e : graph.edges()) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            if (in_handle[u] != in_handle[v]) {
                lhs += ctx.x_values[e];
            }
        }

        // x(delta(T_i)): tooth edges
        for (const auto& tooth : selected) {
            lhs += ctx.x_values[tooth.edge];
        }

        double rhs = std::ceil(3.0 * static_cast<double>(k) / 2.0);

        if (lhs < rhs - tol) {
            Cut cut;

            // Mark tooth edges to avoid duplicates
            std::vector<bool> is_tooth(graph.num_edges(), false);
            for (const auto& tooth : selected) {
                is_tooth[tooth.edge] = true;
            }

            // Handle boundary edges: coefficient -1, or -2 if also a tooth
            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                if (in_handle[u] != in_handle[v]) {
                    double coeff = is_tooth[e] ? -2.0 : -1.0;
                    cut.indices.push_back(ctx.x_offset + e);
                    cut.values.push_back(coeff);
                    is_tooth[e] = false;  // handled
                }
            }

            // Remaining tooth edges (both endpoints inside or outside handle)
            for (const auto& tooth : selected) {
                if (is_tooth[tooth.edge]) {
                    cut.indices.push_back(ctx.x_offset + tooth.edge);
                    cut.values.push_back(-1.0);
                }
            }

            cut.rhs = -rhs;
            cuts.push_back(std::move(cut));
        }
    }

    return cuts;
}

}  // namespace cptp::sep
