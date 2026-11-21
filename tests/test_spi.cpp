#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "preprocess/edge_elimination.h"
#include "sep/separation_context.h"
#include "sep/spi_separator.h"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace {

struct SupportData {
    cptp::digraph graph;
    std::vector<double> capacity;
    std::unique_ptr<cptp::gomory_hu_tree> tree;
};

SupportData build_support(const cptp::Problem& prob, std::span<const double> x_values,
                          double tol = 1e-6) {
    const auto& g = prob.graph();
    cptp::digraph_builder builder(prob.num_nodes());
    for (auto e : g.edges()) {
        const double xval = x_values[e];
        if (xval <= tol) {
            continue;
        }
        const int32_t u = g.edge_source(e);
        const int32_t v = g.edge_target(e);
        builder.add_arc(u, v, xval);
        builder.add_arc(v, u, xval);
    }
    auto [support_graph, cap] = builder.build();
    auto tree = std::make_unique<cptp::gomory_hu_tree>(support_graph, cap, prob.source());
    return {std::move(support_graph), std::move(cap), std::move(tree)};
}

std::vector<double> compute_all_pairs(const cptp::Problem& prob) {
    const int32_t n = prob.num_nodes();
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::vector<double> all_pairs(static_cast<size_t>(n) * n, kInf);
    for (int32_t s = 0; s < n; ++s) {
        auto row = cptp::preprocess::forward_labeling(prob, s);
        std::ranges::copy(row, all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
    }
    return all_pairs;
}

cptp::Problem make_spi_path_problem() {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;

    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            edges.push_back({i, j});
            // Chain edges cheap, shortcuts expensive.
            costs.push_back((j == i + 1) ? 5.0 : 50.0);
        }
    }

    std::vector<double> profits = {0, 0, 0, 0, 0};
    std::vector<double> demands = {0, 1, 1, 1, 0};
    prob.build(5, edges, costs, profits, demands, 10.0, 0, 4);
    return prob;
}

}  // namespace

TEST_CASE("SPI separator: path problem emits violated cuts with all-pairs", "[spi][path]") {
    auto prob = make_spi_path_problem();
    REQUIRE_FALSE(prob.is_tour());

    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    auto all_pairs = compute_all_pairs(prob);

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;  // source
    y_values[4] = 1.0;  // target
    for (int i = 1; i < 4; ++i) {
        y_values[i] = 0.7;
    }

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        const int32_t u = graph.edge_source(e);
        const int32_t v = graph.edge_target(e);
        if (u == 0 || v == 4) {
            x_values[e] = 0.3;
        }
    }

    auto support = build_support(prob, x_values);
    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = 5.0,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    const auto cuts = spi.separate(ctx);
    REQUIRE(!cuts.empty());
    for (const auto& cut : cuts) {
        REQUIRE(cut.violation > 1e-6);
        REQUIRE(cut.rhs >= 1.0 - 1e-9);
        for (double coeff : cut.values) {
            REQUIRE(coeff == 1.0);
        }
    }
}

TEST_CASE("SPI separator: path problem returns no cuts without all-pairs", "[spi][path]") {
    auto prob = make_spi_path_problem();
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[4] = 1.0;
    for (int i = 1; i < 4; ++i) {
        y_values[i] = 0.7;
    }

    auto support = build_support(prob, x_values);
    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = 5.0,
    };

    cptp::sep::SPISeparator spi;
    const auto cuts = spi.separate(ctx);
    REQUIRE(cuts.empty());
}

TEST_CASE("SPI separator: pairs of one conflict clique emit one lifted row", "[spi]") {
    // Depot 0 and customers 1, 2, 3.  A single visit costs 20; two cost at
    // least 40 (shortest paths between customers run through the depot), so
    // under UB 30 every pair conflicts and each of the three violated pairs
    // lifts into y1 + y2 + y3 <= 1 (#176).
    cptp::Problem prob;
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            edges.push_back({i, j});
            costs.push_back(i == 0 ? 10.0 : 100.0);
        }
    }
    const std::vector<double> profits = {0, 0, 0, 0};
    const std::vector<double> demands = {0, 1, 1, 1};
    prob.build(4, edges, costs, profits, demands, 100.0, 0, 0);
    REQUIRE(prob.is_tour());

    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();
    auto all_pairs = compute_all_pairs(prob);

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.6);
    y_values[0] = 1.0;
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        if (graph.edge_source(e) == 0 || graph.edge_target(e) == 0) {
            x_values[e] = 0.6;
        }
    }

    auto support = build_support(prob, x_values);
    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = 30.0,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    const auto cuts = spi.separate(ctx);
    REQUIRE(cuts.size() == 1);
    auto support_cols = cuts[0].indices;
    std::ranges::sort(support_cols);
    REQUIRE(support_cols == std::vector<int32_t>{m + 1, m + 2, m + 3});
    REQUIRE(cuts[0].rhs == 1.0);
}
