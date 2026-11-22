#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "parallel/parallel.h"
#include "sep/largest_side_cache.h"
#include "sep/multistar_separator.h"
#include "sep/separation_context.h"
#include "sep/xcover_separator.h"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

namespace {

/// A complete graph on n nodes with the given demand rows and capacities;
/// source 0, target `target` (0 for a tour).
cptp::Problem make_problem(int32_t n, std::vector<std::vector<double>> rows,
                           std::vector<double> capacities, int32_t target = 0) {
    std::vector<cptp::Edge> edges;
    for (int32_t i = 0; i < n; ++i) {
        for (int32_t j = i + 1; j < n; ++j) {
            edges.push_back({i, j});
        }
    }
    std::vector<double> costs(edges.size(), 1.0);
    std::vector<double> profits(n, 1.0);
    profits[0] = 0.0;
    cptp::Problem prob;
    prob.build(n, edges, std::move(costs), std::move(profits), std::move(rows),
               std::move(capacities), 0, target);
    return prob;
}

bool feasible(const cptp::Problem& prob, const std::vector<int32_t>& nodes) {
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        double load = 0.0;
        for (const int32_t i : nodes) {
            load += prob.demand(k, i);
        }
        if (load > prob.capacity(k) + 1e-9) {
            return false;
        }
    }
    return true;
}

/// Every subset of {1, ..., n-1} as a node list.
std::vector<std::vector<int32_t>> customer_subsets(int32_t n) {
    std::vector<std::vector<int32_t>> subsets;
    for (uint32_t mask = 0; mask < (1U << (n - 1)); ++mask) {
        std::vector<int32_t> nodes;
        for (int32_t i = 1; i < n; ++i) {
            if ((mask >> (i - 1)) & 1U) {
                nodes.push_back(i);
            }
        }
        subsets.push_back(std::move(nodes));
    }
    return subsets;
}

struct Support {
    cptp::digraph graph;
    std::vector<double> capacity;
    std::unique_ptr<cptp::gomory_hu_tree> tree;
};

Support build_support(const cptp::Problem& prob, const std::vector<double>& x) {
    const auto& g = prob.graph();
    cptp::digraph_builder builder(prob.num_nodes());
    for (auto e : g.edges()) {
        if (x[e] > 1e-6) {
            builder.add_arc(g.edge_source(e), g.edge_target(e), x[e]);
            builder.add_arc(g.edge_target(e), g.edge_source(e), x[e]);
        }
    }
    auto [sg, cap] = builder.build();
    auto tree = std::make_unique<cptp::gomory_hu_tree>(sg, cap, prob.source());
    return {std::move(sg), std::move(cap), std::move(tree)};
}

int32_t edge_index(const cptp::Problem& prob, int32_t a, int32_t b) {
    const auto& g = prob.graph();
    for (auto e : g.edges()) {
        const int32_t u = g.edge_source(e);
        const int32_t v = g.edge_target(e);
        if ((u == a && v == b) || (u == b && v == a)) {
            return e;
        }
    }
    return -1;
}

/// (x, y) of every elementary tour through 0 (or s-t path from 0 to the
/// target) that serves a feasible set.
std::vector<std::pair<std::vector<double>, std::vector<double>>> all_routes(
    const cptp::Problem& prob) {
    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();
    const int32_t target = prob.target();
    std::vector<std::pair<std::vector<double>, std::vector<double>>> routes;
    for (auto inner : customer_subsets(n)) {
        if (!prob.is_tour()) {
            if (std::ranges::find(inner, target) == inner.end()) {
                continue;
            }
            std::erase(inner, target);
        }
        std::vector<int32_t> served = inner;
        served.push_back(0);
        if (!prob.is_tour()) {
            served.push_back(target);
        }
        if (!feasible(prob, served)) {
            continue;
        }
        std::ranges::sort(inner);
        do {
            std::vector<int32_t> walk = {0};
            walk.insert(walk.end(), inner.begin(), inner.end());
            walk.push_back(prob.is_tour() ? 0 : target);
            if (prob.is_tour() && inner.empty()) {
                continue;  // the empty tour: x = 0, y_0 = 1
            }
            std::vector<double> x(m, 0.0);
            std::vector<double> y(n, 0.0);
            for (const int32_t v : walk) {
                y[v] = 1.0;
            }
            for (size_t s = 0; s + 1 < walk.size(); ++s) {
                x[edge_index(prob, walk[s], walk[s + 1])] += 1.0;
            }
            routes.emplace_back(std::move(x), std::move(y));
        } while (std::ranges::next_permutation(inner).found);
    }
    return routes;
}

}  // namespace

TEST_CASE("XCover lifting: the #146 example lifts y3, not y4", "[xcover]") {
    // d1 = (1,1,0,1), d2 = (1,1,1,0), Q = (1,1), C = {1,2}, lifting 3 then 4:
    // min over rows gives y1 + y2 + y3 <= 1.
    const auto prob = make_problem(5, {{0, 1, 1, 0, 1}, {0, 1, 1, 1, 0}}, {1.0, 1.0});
    const std::vector<int32_t> cover = {1, 2};
    const std::vector<int32_t> order = {3, 4};
    REQUIRE(cptp::sep::lift_cover_all_rows(prob, cover, order) ==
            std::vector<int32_t>{0, 1, 1, 1, 0});
}

TEST_CASE("XCover lifting: every feasible set keeps the lifted cover", "[xcover]") {
    std::mt19937 rng(166);
    std::uniform_int_distribution<int32_t> demand(0, 5);
    std::uniform_int_distribution<int32_t> capacity(5, 10);
    constexpr int32_t kNodes = 8;
    const auto subsets = customer_subsets(kNodes);
    int32_t checked = 0;
    for (int32_t trial = 0; trial < 60; ++trial) {
        const int32_t m = 1 + (trial % 3);
        std::vector<std::vector<double>> rows(m, std::vector<double>(kNodes, 0.0));
        std::vector<double> caps(m);
        for (int32_t k = 0; k < m; ++k) {
            caps[k] = capacity(rng);
            for (int32_t i = 1; i < kNodes; ++i) {
                rows[k][i] = demand(rng);
            }
        }
        const auto prob = make_problem(kNodes, rows, caps);
        for (const auto& cover : subsets) {
            if (cover.size() < 2) {
                continue;
            }
            bool is_cover = false;
            for (int32_t k = 0; k < m; ++k) {
                double load = 0.0;
                for (const int32_t i : cover) {
                    load += rows[k][i];
                }
                is_cover = is_cover || load > caps[k];
            }
            if (!is_cover) {
                continue;
            }
            std::vector<int32_t> order;
            for (int32_t i = 1; i < kNodes; ++i) {
                if (std::ranges::find(cover, i) == cover.end()) {
                    order.push_back(i);
                }
            }
            std::ranges::shuffle(order, rng);
            const auto alpha = cptp::sep::lift_cover_all_rows(prob, cover, order);
            const auto rhs = static_cast<int32_t>(cover.size()) - 1;
            for (const auto& set : subsets) {
                if (!feasible(prob, set)) {
                    continue;
                }
                int32_t lhs = 0;
                for (const int32_t i : set) {
                    lhs += alpha[i];
                }
                INFO("trial " << trial);
                REQUIRE(lhs <= rhs);
            }
            ++checked;
        }
    }
    REQUIRE(checked > 100);
}

TEST_CASE("XCover cuts: no feasible tour or s-t path violates them", "[xcover]") {
    std::mt19937 rng(1660);
    std::uniform_int_distribution<int32_t> demand(1, 4);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    constexpr int32_t kNodes = 6;
    int32_t cuts_checked = 0;
    for (int32_t trial = 0; trial < 40; ++trial) {
        const int32_t target = trial % 2 == 0 ? 0 : kNodes - 1;
        std::vector<std::vector<double>> rows(2, std::vector<double>(kNodes, 0.0));
        for (auto& row : rows) {
            for (int32_t i = 1; i < kNodes; ++i) {
                row[i] = demand(rng);
            }
        }
        const auto prob = make_problem(kNodes, rows, {6.0, 6.0}, target);
        const auto routes = all_routes(prob);
        REQUIRE_FALSE(routes.empty());
        const int32_t m = prob.num_edges();

        // A random fractional point, every node partly served.
        std::vector<double> x(m);
        std::vector<double> y(kNodes);
        for (auto& v : x) {
            v = unit(rng) < 0.5 ? 0.0 : unit(rng);
        }
        for (auto& v : y) {
            v = 0.3 + 0.7 * unit(rng);
        }
        y[0] = 1.0;
        y[target] = 1.0;
        auto support = build_support(prob, x);
        const cptp::sep::SeparationContext ctx{
            .problem = prob,
            .x_values = x,
            .y_values = y,
            .x_offset = 0,
            .y_offset = m,
            .tol = 1e-6,
            .flow_tree = support.tree.get(),
            .support_graph = &support.graph,
            .support_capacity = &support.capacity,
        };
        cptp::sep::XCoverSeparator xcover;
        for (const auto& cut : xcover.separate(ctx)) {
            REQUIRE(cut.violation > 1e-6);
            for (const auto& [rx, ry] : routes) {
                double lhs = 0.0;
                for (int32_t t = 0; t < cut.size(); ++t) {
                    const int32_t col = cut.indices[t];
                    lhs += cut.values[t] * (col < m ? rx[col] : ry[col - m]);
                }
                INFO("trial " << trial);
                REQUIRE(lhs <= cut.rhs + 1e-9);
            }
            ++cuts_checked;
        }
    }
    REQUIRE(cuts_checked > 20);
}

TEST_CASE("Largest sides: one max-flow per target for multistar and covers together",
          "[xcover][multistar]") {
    // Both separators ask for every target with y > tol; run in parallel on a
    // shared cache they cost one depot-target max-flow per target, not two,
    // and give the cuts each gives on its own max-flows (#175).
    std::mt19937 rng(175);
    std::uniform_int_distribution<int32_t> demand(1, 4);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    constexpr int32_t kNodes = 8;
    const auto same = [](const std::vector<cptp::sep::Cut>& a,
                         const std::vector<cptp::sep::Cut>& b) {
        return std::ranges::equal(a, b, [](const cptp::sep::Cut& p, const cptp::sep::Cut& q) {
            return p.indices == q.indices && p.values == q.values && p.rhs == q.rhs &&
                   p.violation == q.violation;
        });
    };
    int32_t cuts_compared = 0;
    for (int32_t trial = 0; trial < 40; ++trial) {
        const int32_t target = trial % 2 == 0 ? 0 : kNodes - 1;
        std::vector<std::vector<double>> rows(2, std::vector<double>(kNodes, 0.0));
        for (auto& row : rows) {
            for (int32_t i = 1; i < kNodes; ++i) {
                row[i] = demand(rng);
            }
        }
        const auto prob = make_problem(kNodes, rows, {7.0, 7.0}, target);
        const int32_t m = prob.num_edges();
        std::vector<double> x(m);
        std::vector<double> y(kNodes);
        // Few distinct values, so that min cuts tie.
        for (auto& v : x) {
            v = unit(rng) < 0.5 ? 0.0 : 0.25 * static_cast<double>(1 + (rng() % 4));
        }
        for (auto& v : y) {
            v = unit(rng) < 0.2 ? 0.0 : 0.5 + 0.5 * unit(rng);
        }
        y[0] = 1.0;
        y[target] = 1.0;
        auto support = build_support(prob, x);
        cptp::sep::SeparationContext ctx{
            .problem = prob,
            .x_values = x,
            .y_values = y,
            .x_offset = 0,
            .y_offset = m,
            .tol = 1e-6,
            .flow_tree = support.tree.get(),
            .support_graph = &support.graph,
            .support_capacity = &support.capacity,
        };
        int64_t targets = 0;
        for (int32_t i = 1; i < kNodes; ++i) {
            targets += (i != target && y[i] > ctx.tol) ? 1 : 0;
        }

        cptp::sep::MultistarSeparator multistar(cptp::sep::CapacityRows::kAll, true);
        cptp::sep::XCoverSeparator xcover;
        const auto multistar_alone = multistar.separate(ctx);
        const auto xcover_alone = xcover.separate(ctx);

        cptp::sep::LargestSideCache shared(support.graph, support.capacity, prob.source());
        ctx.largest_sides = &shared;
        std::vector<cptp::sep::Cut> multistar_shared;
        std::vector<cptp::sep::Cut> xcover_shared;
        {
            cptp::parallel::task_group tg;
            tg.run([&] { multistar_shared = multistar.separate(ctx); });
            tg.run([&] { xcover_shared = xcover.separate(ctx); });
        }
        INFO("trial " << trial);
        REQUIRE(shared.max_flows() == targets);
        REQUIRE(same(multistar_shared, multistar_alone));
        REQUIRE(same(xcover_shared, xcover_alone));
        cuts_compared += static_cast<int32_t>(multistar_alone.size() + xcover_alone.size());

        // A second round on the same cache runs no max-flow.
        REQUIRE(same(multistar.separate(ctx), multistar_alone));
        REQUIRE(shared.max_flows() == targets);
    }
    REQUIRE(cuts_compared > 20);
}
