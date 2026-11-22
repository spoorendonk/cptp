#include "sep/xcover_separator.h"

#include "core/dinitz.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "sep/largest_side_cache.h"

#include <algorithm>
#include <limits>
#include <optional>
#include <set>
#include <vector>

namespace cptp::sep {

namespace {

/// Slack on capacity comparisons, relative to the capacity.
double capacity_eps(double Q) {
    return 1e-9 * std::max(1.0, Q);
}

/// A cover of row `row` over the candidates, greedy by (1 - z_i) / d_i as in
/// Crowder, Johnson and Padberg's separation heuristic, then made minimal by
/// dropping the smallest z first. Empty when the candidates fit.
std::vector<int32_t> greedy_cover(const Problem& prob, int32_t row,
                                  const std::vector<int32_t>& candidates,
                                  const std::vector<double>& z) {
    const double Q = prob.capacity(row);
    const double limit = Q + capacity_eps(Q);
    std::vector<int32_t> items;
    for (const int32_t i : candidates) {
        if (prob.demand(row, i) > 0) {
            items.push_back(i);
        }
    }
    std::ranges::sort(items, [&](int32_t a, int32_t b) {
        const double ra = (1.0 - z[a]) / prob.demand(row, a);
        const double rb = (1.0 - z[b]) / prob.demand(row, b);
        return ra != rb ? ra < rb : a < b;
    });
    std::vector<int32_t> cover;
    double load = 0.0;
    for (const int32_t i : items) {
        cover.push_back(i);
        load += prob.demand(row, i);
        if (load > limit) {
            break;
        }
    }
    if (load <= limit) {
        return {};
    }
    std::ranges::sort(cover,
                      [&](int32_t a, int32_t b) { return z[a] != z[b] ? z[a] < z[b] : a < b; });
    std::vector<int32_t> minimal;
    for (const int32_t i : cover) {
        if (load - prob.demand(row, i) > limit) {
            load -= prob.demand(row, i);
        } else {
            minimal.push_back(i);
        }
    }
    std::ranges::sort(minimal);
    return minimal;
}

}  // namespace

std::vector<int32_t> lift_cover_all_rows(const Problem& prob, std::span<const int32_t> cover,
                                         std::span<const int32_t> order) {
    const int32_t n = prob.num_nodes();
    const auto rhs = static_cast<int32_t>(cover.size()) - 1;
    std::vector<int32_t> alpha(n, 0);
    std::vector<int32_t> items;
    for (const int32_t i : cover) {
        alpha[i] = 1;
        items.push_back(i);
    }
    std::vector<int32_t> rows;
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        if (prob.capacity(k) < kInfiniteCapacity) {
            rows.push_back(k);
        }
    }

    // min_weight[v]: the least row-k demand of a subset of `items` whose
    // value, capped at rhs, is v. The exact lifting value never exceeds rhs,
    // so the cap loses nothing.
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::vector<double> min_weight(rhs + 1);
    for (const int32_t j : order) {
        int32_t bound = rhs;  // U_j: least single-row lifting value so far
        bool alone_infeasible = false;
        for (const int32_t k : rows) {
            const double Q = prob.capacity(k);
            const double room = Q - prob.demand(k, j);
            if (room < -capacity_eps(Q)) {
                alone_infeasible = true;
                break;
            }
            std::ranges::fill(min_weight, kInf);
            min_weight[0] = 0.0;
            for (const int32_t i : items) {
                const int32_t a = alpha[i];
                const double d = prob.demand(k, i);
                for (int32_t v = rhs; v >= 0; --v) {
                    if (min_weight[v] == kInf) {
                        continue;
                    }
                    const int32_t nv = std::min(rhs, v + a);
                    min_weight[nv] = std::min(min_weight[nv], min_weight[v] + d);
                }
            }
            int32_t best = 0;
            for (int32_t v = rhs; v > 0; --v) {
                if (min_weight[v] <= room + capacity_eps(Q)) {
                    best = v;
                    break;
                }
            }
            bound = std::min(bound, best);
        }
        // No feasible set holds a node that cannot be served alone, so it
        // stays out of the later lifting problems.
        alpha[j] = alone_infeasible ? rhs : rhs - bound;
        if (!alone_infeasible && alpha[j] > 0) {
            items.push_back(j);
        }
    }
    return alpha;
}

std::vector<Cut> XCoverSeparator::separate(const SeparationContext& ctx) {
    const auto& prob = ctx.problem;
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const double tol = ctx.tol;
    const int32_t depot = prob.source();
    const bool is_tour = prob.is_tour();
    const int32_t path_target = prob.target();

    std::vector<int32_t> cover_rows;
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        const double Q = prob.capacity(k);
        if (Q > 0 && Q < kInfiniteCapacity) {
            cover_rows.push_back(k);
        }
    }
    if (cover_rows.empty()) {
        return {};
    }

    // Candidate sets S, as membership vectors without the depot and, for an
    // s-t path, without its target (as the multistar separator's GLMs).
    std::set<std::vector<bool>> sets;
    auto add_set = [&](const std::vector<bool>& side) {
        std::vector<bool> in_S(n);
        bool any = false;
        for (int32_t u = 0; u < n; ++u) {
            in_S[u] = side[u] && u != depot && (is_tour || u != path_target);
            any = any || in_S[u];
        }
        if (any) {
            sets.insert(std::move(in_S));
        }
    };
    const bool have_support = ctx.support_graph != nullptr && ctx.support_capacity != nullptr;
    std::optional<LargestSideCache> own_sides;
    LargestSideCache* largest_sides = ctx.largest_sides;
    if (have_support && largest_sides == nullptr) {
        largest_sides = &own_sides.emplace(*ctx.support_graph, *ctx.support_capacity, depot);
    }
    for (int32_t target = 0; target < n; ++target) {
        if (target == depot || (!is_tour && target == path_target) || ctx.y_values[target] <= tol) {
            continue;
        }
        const auto reachable = ctx.flow_tree->min_cut(target).on_source_side;
        std::vector<bool> smallest(n);
        for (int32_t u = 0; u < n; ++u) {
            smallest[u] = !reachable[u];
        }
        add_set(smallest);
        if (have_support) {
            add_set(largest_sides->get(target).on_sink_side);
        }
    }
    add_set(std::vector<bool>(n, true));

    // A boundary edge of S and its endpoint outside S, found once per S for
    // z and every row's cover (#175).
    struct Crossing {
        int32_t edge;
        int32_t outside;
    };
    std::vector<Crossing> boundary;

    std::vector<Cut> cuts;
    std::vector<double> z(n);
    for (const auto& in_S : sets) {
        boundary.clear();
        for (auto e : graph.edges()) {
            const int32_t u = graph.edge_source(e);
            const int32_t v = graph.edge_target(e);
            if (in_S[u] != in_S[v]) {
                boundary.push_back({e, in_S[u] ? v : u});
            }
        }

        // z_i = y_i in S, z_j = x(j : S) outside (the depot left out).
        std::ranges::fill(z, 0.0);
        for (int32_t i = 0; i < n; ++i) {
            if (in_S[i]) {
                z[i] = ctx.y_values[i];
            }
        }
        for (const auto [e, outside] : boundary) {
            if (outside != depot) {
                z[outside] += ctx.x_values[e];
            }
        }
        std::vector<int32_t> candidates;
        for (int32_t i = 0; i < n; ++i) {
            if (i != depot && z[i] > 1e-9) {
                candidates.push_back(i);
            }
        }

        std::set<std::vector<int32_t>> covers_seen;
        for (const int32_t row : cover_rows) {
            auto cover = greedy_cover(prob, row, candidates, z);
            if (cover.size() < 2 || !covers_seen.insert(cover).second) {
                continue;
            }
            std::vector<int32_t> order;
            for (const int32_t i : candidates) {
                if (!std::ranges::binary_search(cover, i)) {
                    order.push_back(i);
                }
            }
            std::ranges::sort(
                order, [&](int32_t a, int32_t b) { return z[a] != z[b] ? z[a] > z[b] : a < b; });
            const auto alpha = lift_cover_all_rows(prob, cover, order);
            const auto rhs = static_cast<double>(cover.size() - 1);

            // w-GLM in <= form, w = alpha / rhs:
            //   2 sum_{i in S} w_i y_i - sum_{e in delta(S)} (1 - 2 w_{t(e)}) x_e <= 0
            Cut cut;
            double lhs = 0.0;
            for (const auto [e, outside] : boundary) {
                const double w = outside == depot ? 0.0 : alpha[outside] / rhs;
                const double coeff = 1.0 - 2.0 * w;
                lhs -= coeff * ctx.x_values[e];
                cut.indices.push_back(ctx.x_offset + e);
                cut.values.push_back(-coeff);
            }
            for (int32_t i = 0; i < n; ++i) {
                if (!in_S[i] || alpha[i] == 0) {
                    continue;
                }
                const double coeff = 2.0 * alpha[i] / rhs;
                lhs += coeff * ctx.y_values[i];
                cut.indices.push_back(ctx.y_offset + i);
                cut.values.push_back(coeff);
            }
            cut.rhs = 0.0;
            if (lhs <= tol || cut.indices.empty()) {
                continue;
            }
            cut.violation = lhs;
            cuts.push_back(std::move(cut));
        }
    }
    return cuts;
}

}  // namespace cptp::sep
