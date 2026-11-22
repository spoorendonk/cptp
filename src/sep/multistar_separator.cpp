#include "sep/multistar_separator.h"

#include "core/dinitz.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "sep/capacity_rows.h"
#include "sep/largest_side_cache.h"

#include <cmath>
#include <map>
#include <optional>
#include <utility>
#include <vector>

namespace cptp::sep {

std::vector<Cut> MultistarSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.source();
    const bool is_tour = prob.is_tour();
    const int32_t path_target = prob.target();

    // Rows that can give a cut: 0 < Q_k < kInfiniteCapacity, as row 0 alone
    // was before, with capacity left after the terminals' demand (#170).
    struct Row {
        int32_t index;
        double residual;
    };
    std::vector<Row> rows;
    const int32_t num_rows = rows_ == CapacityRows::kFirst ? 1 : prob.num_capacities();
    for (int32_t k = 0; k < num_rows; ++k) {
        const double Q = prob.capacity(k);
        if (Q <= 0 || Q >= kInfiniteCapacity) {
            continue;
        }
        const double residual = residual_capacity(prob, {prob.demands(k), Q});
        if (residual > 0.0) {
            rows.push_back({k, residual});
        }
    }
    if (rows.empty()) {
        return {};
    }

    std::vector<Cut> all_cuts;

    // A boundary edge of S and its endpoint outside S.
    struct Crossing {
        int32_t edge;
        int32_t outside;
    };
    std::vector<Crossing> boundary;

    // The GLMs on S (the side without the depot, and without the path's
    // target), per row k (d = d^k, Q = Q_k less the terminals' demand, which
    // is weighted 0: see residual_capacity):
    //   sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e
    //     - (2/Q)*sum_{i in S} q_i*y_i >= 0
    // In <= form:
    //   (2/Q)*sum_{i in S} q_i*y_i
    //     - sum_{e in delta(S)} (1 - 2*q_{t(e)}/Q)*x_e <= 0
    // One cut per row (kAll), or the most violated row's (otherwise).
    // delta(S) is the same for every row, so it is found once (#175).
    auto add_glms = [&](const std::vector<bool>& in_S) {
        boundary.clear();
        for (auto e : graph.edges()) {
            const int32_t u = graph.edge_source(e);
            const int32_t v = graph.edge_target(e);
            if (in_S[u] != in_S[v]) {
                boundary.push_back({e, in_S[u] ? v : u});
            }
        }

        Cut best;
        bool have_best = false;

        for (const auto [row, Q] : rows) {
            double lhs = 0.0;
            Cut cut;

            for (const auto [e, outside] : boundary) {
                const double q_outside =
                    is_terminal(prob, outside) ? 0.0 : prob.demand(row, outside);
                const double coeff = 1.0 - 2.0 * q_outside / Q;
                lhs -= coeff * ctx.x_values[e];
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-coeff);
            }

            for (int32_t i = 0; i < n; ++i) {
                if (!in_S[i]) {
                    continue;
                }
                const double q_i = prob.demand(row, i);
                if (q_i <= 0) {
                    continue;
                }

                const double coeff = 2.0 * q_i / Q;
                lhs += coeff * ctx.y_values[i];
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(coeff);
            }

            cut.rhs = 0.0;

            if (lhs <= tol || cut.indices.empty()) {
                continue;
            }
            cut.violation = lhs;
            if (rows_ == CapacityRows::kAll) {
                all_cuts.push_back(std::move(cut));
            } else if (!have_best || cut.violation > best.violation) {
                best = std::move(cut);
                have_best = true;
            }
        }

        if (have_best) {
            all_cuts.push_back(std::move(best));
        }
    };

    // A side as the set S its GLMs are on: the path target stripped (the
    // depot is never on a target's side).
    auto effective = [&](const std::vector<bool>& side) {
        std::vector<bool> in_S(side);
        if (!is_tour) {
            in_S[path_target] = false;
        }
        return in_S;
    };

    // The sets S cut so far, each with its cuts' range in all_cuts. Targets
    // sharing a Gomory-Hu bottleneck give the same S. With largest_side each
    // S is cut once over all targets, smallest or largest side. Without it
    // the smallest sides are cut as before, even when repeated: a repeat
    // copies the first one's cuts rather than building them again (#175).
    std::map<std::vector<bool>, std::pair<size_t, size_t>> cut_sets;
    auto cut_set = [&](std::vector<bool> in_S, bool repeat) {
        auto [it, inserted] = cut_sets.try_emplace(std::move(in_S));
        if (!inserted) {
            if (repeat) {
                for (size_t c = it->second.first; c < it->second.second; ++c) {
                    Cut copy = all_cuts[c];
                    all_cuts.push_back(std::move(copy));
                }
            }
            return;
        }
        const size_t begin = all_cuts.size();
        add_glms(it->first);
        it->second = {begin, all_cuts.size()};
    };

    const bool try_largest =
        largest_side_ && ctx.support_graph != nullptr && ctx.support_capacity != nullptr;
    std::optional<LargestSideCache> own_sides;
    LargestSideCache* largest_sides = ctx.largest_sides;
    if (try_largest && largest_sides == nullptr) {
        largest_sides = &own_sides.emplace(*ctx.support_graph, *ctx.support_capacity, depot);
    }

    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || (!is_tour && target == path_target) || ctx.y_values[target] <= tol) {
            continue;
        }

        auto [flow, reachable] = ctx.flow_tree->min_cut(target);
        std::vector<bool> smallest(n);
        for (int32_t u = 0; u < n; ++u) {
            smallest[u] = !reachable[u];
        }
        cut_set(effective(smallest), !try_largest);
        if (try_largest) {
            cut_set(effective(largest_sides->get(target).on_sink_side), false);
        }
    }

    return all_cuts;
}

}  // namespace cptp::sep
