#include "sep/rci_separator.h"

#include <cmath>

#include "core/digraph.h"
#include "core/dinitz.h"
#include "core/problem.h"

#include <tbb/task_group.h>
#include <mutex>

namespace cptp::sep {

std::vector<Cut> RCISeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const double Q = prob.capacity();
    const int32_t depot = prob.depot();

    if (Q <= 0) return {};

    // Build directed support graph
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
            double flow = alg.flow_value();

            // Get the min-cut set S = nodes NOT reachable from depot
            std::vector<bool> reachable(n, false);
            for (int32_t u = 0; u < n; ++u) {
                reachable[u] = alg.on_source_side(u);
            }

            // Compute total demand of S (unweighted)
            double d_S = 0.0;
            for (int32_t i = 0; i < n; ++i) {
                if (!reachable[i]) {
                    d_S += prob.demand(i);
                }
            }

            double Q_r = std::fmod(d_S, Q);
            if (Q_r <= tol) return;  // d(S) is multiple of Q, no strengthening

            double k = std::ceil(d_S / Q);
            // Only useful when k > 1 (otherwise same as SEC)
            if (k <= 1.0) return;

            // CPTP RCI (Jepsen et al. 2014):
            //   x(delta(S)) - (2/Q_r) * sum_{i in S} q_i * y_i
            //     >= 2*(ceil(d(S)/Q) - d(S)/Q_r)
            //
            // Compute weighted demand served: sum_{i in S} q_i * y_i
            double served_demand = 0.0;
            for (int32_t i = 0; i < n; ++i) {
                if (!reachable[i]) {
                    served_demand += prob.demand(i) * ctx.y_values[i];
                }
            }

            double rhs = 2.0 * (k - d_S / Q_r);
            double lhs = flow - (2.0 * served_demand) / Q_r;

            if (lhs >= rhs - tol) return;  // not violated

            // Build the cut in <= form:
            //   (2/Q_r) * sum_{i in S} q_i * y_i - x(delta(S)) <= -rhs
            //   i.e., -rhs is our cut.rhs
            Cut cut;

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
                std::lock_guard lock(cuts_mutex);
                all_cuts.push_back(std::move(cut));
            }
        });
    }
    tg.wait();

    return all_cuts;
}

}  // namespace cptp::sep
