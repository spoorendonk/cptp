#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/io.h"
#include "core/problem.h"
#include "sep/sec_separator.h"
#include "sep/rci_separator.h"
#include "sep/separation_context.h"

using Catch::Matchers::WithinAbs;

namespace {

/// Build a small 4-node problem: depot=0, customers=1,2,3.
/// Complete undirected graph (6 edges for 4 nodes).
cptp::Problem make_small_problem() {
    cptp::Problem prob;

    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            edges.push_back({i, j});
            costs.push_back(10.0);  // uniform cost
        }
    }
    // Edges: {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}

    std::vector<double> profits = {0, 20, 15, 10};
    std::vector<double> demands = {0, 3, 4, 2};
    double capacity = 7.0;

    prob.build(4, edges, costs, profits, demands, capacity, 0, 0);
    return prob;
}

/// Build a small 4-node s-t path problem: source=0, target=3.
cptp::Problem make_small_path_problem() {
    cptp::Problem prob;

    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            edges.push_back({i, j});
            costs.push_back(10.0);
        }
    }

    std::vector<double> profits = {0, 20, 15, 10};
    std::vector<double> demands = {0, 3, 4, 2};
    double capacity = 7.0;

    prob.build(4, edges, costs, profits, demands, capacity, 0, 3);
    return prob;
}

/// Build support graph and Gomory-Hu tree from LP solution.
struct SupportData {
    cptp::digraph graph;
    std::vector<double> capacity;
    std::unique_ptr<cptp::gomory_hu_tree> tree;
};

SupportData build_support(const cptp::Problem& prob,
                          std::span<const double> x_values,
                          double tol = 1e-6) {
    const auto& g = prob.graph();
    int32_t n = prob.num_nodes();
    cptp::digraph_builder builder(n);
    for (auto e : g.edges()) {
        double xval = x_values[e];
        if (xval > tol) {
            int32_t u = g.edge_source(e);
            int32_t v = g.edge_target(e);
            builder.add_arc(u, v, xval);
            builder.add_arc(v, u, xval);
        }
    }
    auto [sg, cap] = builder.build();
    auto tree = std::make_unique<cptp::gomory_hu_tree>(sg, cap, prob.source());
    return {std::move(sg), std::move(cap), std::move(tree)};
}

}  // namespace

TEST_CASE("SEC separator finds violated cuts on fractional solution", "[sec]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Create a fractional LP solution that violates SEC.
    // Subtour {0,1} and disconnected subtour {2,3}
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;  // depot
    y_values[1] = 0.5;
    y_values[2] = 0.5;
    y_values[3] = 0.5;

    // Edges: {0,1}=0, {0,2}=1, {0,3}=2, {1,2}=3, {1,3}=4, {2,3}=5
    // Set fractional edges: {0,1} at 1.0 (degree constraint: 2*y_i)
    // {2,3} at 1.0 (disconnected subtour)
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 1) {
            x_values[e] = 1.0;
        }
        if (u == 2 && v == 3) {
            x_values[e] = 1.0;
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
    };

    cptp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // Should find violated SEC for nodes 2 and/or 3 (disconnected from depot)
    REQUIRE(!cuts.empty());
}

TEST_CASE("SEC separator finds no cuts on integer feasible solution", "[sec]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    // Feasible tour: 0 - 1 - 2 - 0, node 3 not visited
    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[2] = 1.0;
    y_values[3] = 0.0;

    // Edges: {0,1}=0, {0,2}=1, {1,2}=3
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 0 && v == 2) || (u == 1 && v == 2)) {
            x_values[e] = 1.0;
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
    };

    cptp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // No violated cuts on a feasible solution
    REQUIRE(cuts.empty());
}

TEST_CASE("RCI separator basic test", "[rci]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[2] = 1.0;
    y_values[3] = 0.0;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 0 && v == 2) || (u == 1 && v == 2)) {
            x_values[e] = 1.0;
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
    };

    cptp::sep::RCISeparator rci;
    auto cuts = rci.separate(ctx);

    // Just verify it doesn't crash
    REQUIRE(true);
}

TEST_CASE("SEC separator for s-t path finds violated cuts", "[sec][path]") {
    auto prob = make_small_path_problem();
    REQUIRE_FALSE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 3);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Create a solution where node 2 is disconnected from source
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;  // source
    y_values[3] = 1.0;  // target
    y_values[1] = 0.5;
    y_values[2] = 0.5;

    // Only connect 0-3 and 1-2 (node 2 disconnected from source)
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 3) || (u == 1 && v == 2)) {
            x_values[e] = 1.0;
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
    };

    cptp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // Should find violated SEC for disconnected nodes
    REQUIRE(!cuts.empty());
}

TEST_CASE("SEC separator for s-t path: feasible path has no violations", "[sec][path]") {
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Feasible path: 0 -> 1 -> 3 (source=0, target=3)
    // Degree: source=1, node1=2, target=1
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;  // source
    y_values[1] = 1.0;  // intermediate
    y_values[3] = 1.0;  // target

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 1 && v == 3)) {
            x_values[e] = 1.0;
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
    };

    cptp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // No violated cuts on a feasible path
    REQUIRE(cuts.empty());
}

TEST_CASE("Problem is_tour() and source/target accessors", "[problem]") {
    cptp::Problem tour_prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 1, 1};

    tour_prob.build(3, edges, costs, profits, demands, 10.0, 0, 0);
    REQUIRE(tour_prob.is_tour());
    REQUIRE(tour_prob.source() == 0);
    REQUIRE(tour_prob.target() == 0);
    REQUIRE(tour_prob.depot() == 0);

    cptp::Problem path_prob;
    path_prob.build(3, edges, costs, profits, demands, 10.0, 0, 2);
    REQUIRE_FALSE(path_prob.is_tour());
    REQUIRE(path_prob.source() == 0);
    REQUIRE(path_prob.target() == 2);
    REQUIRE(path_prob.depot() == 0);  // backward compat: returns source
}

TEST_CASE("PathWyse IO: load tour (no source/target line)", "[io]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    REQUIRE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 0);
    REQUIRE(prob.num_nodes() == 4);
    REQUIRE(prob.capacity() == 7.0);
}

TEST_CASE("PathWyse IO: load path (with source/target line)", "[io]") {
    auto prob = cptp::io::load("tests/data/tiny4_path.txt");
    REQUIRE_FALSE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 3);
    REQUIRE(prob.num_nodes() == 4);
    REQUIRE(prob.capacity() == 7.0);
}
