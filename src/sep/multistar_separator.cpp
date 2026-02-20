#include "sep/multistar_separator.h"

#include <cmath>

#include "core/digraph.h"
#include "core/dinitz.h"
#include "core/problem.h"

#include <tbb/task_group.h>
#include <mutex>

namespace cptp::sep {

std::vector<Cut> MultistarSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();

    if (Q <= 0 || Q >= 1e17) return {};

    // Build directed support graph (same as SEC)
    digraph_builder builder(n);

    for (auto e : graph.edges()) {
        double xval = ctx.x_values[e];
        if (xval > tol) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            builder.add_arc(u, v, xval);
            builder.add_arc(v, u, xval);
        }
    }

    auto [support_graph, capacity] = builder.build();

    std::vector<Cut> all_cuts;
    std::mutex cuts_mutex;

    // Collect candidate nodes
    std::vector<int32_t> candidates;
    for (int32_t i = 0; i < n; ++i) {
        if (i != depot && ctx.y_values[i] > tol) {
            candidates.push_back(i);
        }
    }

    tbb::task_group tg;
    for (auto target : candidates) {
        tg.run([&, target] {
            dinitz alg(support_graph, capacity, depot, target);
            alg.run();

            // Get the min-cut set S = nodes NOT reachable from depot
            // (i.e., target side of the cut)
            std::vector<bool> reachable(n, false);
            for (int32_t u = 0; u < n; ++u) {
                reachable[u] = alg.on_source_side(u);
            }

            // S = {u : !reachable[u]}, should not contain depot
            // Compute GLM (Generalized Large Multistar) for CPTP:
            //   sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q) * x_e
            //     - (2/Q) * sum_{i in S} q_i * y_i >= 0
            // where t(e) is the endpoint of e NOT in S.
            //
            // In <= form for our Cut struct:
            //   (2/Q) * sum_{i in S} q_i * y_i
            //     - sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q) * x_e <= 0

            double lhs = 0.0;
            Cut cut;

            for (auto e : graph.edges()) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                bool u_in_S = !reachable[u];
                bool v_in_S = !reachable[v];
                if (u_in_S == v_in_S) continue;  // not a crossing edge

                // t(e) = endpoint NOT in S
                int32_t outside = u_in_S ? v : u;
                double q_outside = prob.demand(outside);

                double coeff = 1.0 - 2.0 * q_outside / Q;
                lhs -= coeff * ctx.x_values[e];
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-coeff);
            }

            // y terms for nodes in S
            for (int32_t i = 0; i < n; ++i) {
                if (reachable[i]) continue;  // not in S
                double q_i = prob.demand(i);
                if (q_i <= 0) continue;

                double coeff = 2.0 * q_i / Q;
                lhs += coeff * ctx.y_values[i];
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(coeff);
            }

            cut.rhs = 0.0;

            if (lhs > tol && !cut.indices.empty()) {
                std::lock_guard lock(cuts_mutex);
                all_cuts.push_back(std::move(cut));
            }
        });
    }
    tg.wait();

    return all_cuts;
}

}  // namespace cptp::sep
