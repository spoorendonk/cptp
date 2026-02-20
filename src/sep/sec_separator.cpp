#include "sep/sec_separator.h"

#include "core/digraph.h"
#include "core/dinitz.h"
#include "core/problem.h"

#include <tbb/task_group.h>
#include <mutex>

namespace cptp::sep {

std::vector<Cut> SECSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.depot();

    // Build directed support graph: for each undirected edge {u,v} with x_e > tol,
    // add two arcs (u->v, x_e) and (v->u, x_e)
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

    // For each customer node with y_i > tol, compute min-cut from depot
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
            double yi = ctx.y_values[target];

            // Undirected SEC: flow from depot to i >= 2*y_i
            if (flow < 2.0 * yi - tol) {
                // Determine set S (depot side) from the min-cut.
                std::vector<bool> reachable(n, false);
                for (int32_t u = 0; u < n; ++u) {
                    reachable[u] = alg.on_source_side(u);
                }

                // S = reachable (source side), V\S = not reachable
                // SEC: sum(x_e for edges crossing S to V\S) >= 2*y_target
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
                    std::lock_guard lock(cuts_mutex);
                    all_cuts.push_back(std::move(cut));
                }
            }
        });
    }
    tg.wait();

    return all_cuts;
}

}  // namespace cptp::sep
