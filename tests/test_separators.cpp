#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/io.h"
#include "core/problem.h"
#include "heuristic/warm_start.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/reachability.h"
#include "sep/sec_separator.h"
#include "sep/rci_separator.h"
#include "sep/separation_context.h"

using Catch::Matchers::WithinAbs;

namespace {

/// Build a small 4-node problem: depot=0, customers=1,2,3.
/// Complete undirected graph (6 edges for 4 nodes).
rcspp::Problem make_small_problem() {
    rcspp::Problem prob;

    std::vector<rcspp::Edge> edges;
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
rcspp::Problem make_small_path_problem() {
    rcspp::Problem prob;

    std::vector<rcspp::Edge> edges;
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
    rcspp::digraph graph;
    std::vector<double> capacity;
    std::unique_ptr<rcspp::gomory_hu_tree> tree;
};

SupportData build_support(const rcspp::Problem& prob,
                          std::span<const double> x_values,
                          double tol = 1e-6) {
    const auto& g = prob.graph();
    int32_t n = prob.num_nodes();
    rcspp::digraph_builder builder(n);
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
    auto tree = std::make_unique<rcspp::gomory_hu_tree>(sg, cap, prob.source());
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

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
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

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
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

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::RCISeparator rci;
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

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
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

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // No violated cuts on a feasible path
    REQUIRE(cuts.empty());
}

TEST_CASE("Problem is_tour() and source/target accessors", "[problem]") {
    rcspp::Problem tour_prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 1, 1};

    tour_prob.build(3, edges, costs, profits, demands, 10.0, 0, 0);
    REQUIRE(tour_prob.is_tour());
    REQUIRE(tour_prob.source() == 0);
    REQUIRE(tour_prob.target() == 0);
    REQUIRE(tour_prob.depot() == 0);

    rcspp::Problem path_prob;
    path_prob.build(3, edges, costs, profits, demands, 10.0, 0, 2);
    REQUIRE_FALSE(path_prob.is_tour());
    REQUIRE(path_prob.source() == 0);
    REQUIRE(path_prob.target() == 2);
    REQUIRE(path_prob.depot() == 0);  // backward compat: returns source
}

TEST_CASE("PathWyse IO: load tour (no source/target line)", "[io]") {
    auto prob = rcspp::io::load("tests/data/tiny4.txt");
    REQUIRE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 0);
    REQUIRE(prob.num_nodes() == 4);
    REQUIRE(prob.capacity() == 7.0);
}

TEST_CASE("PathWyse IO: load path (with source/target line)", "[io]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");
    REQUIRE_FALSE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 3);
    REQUIRE(prob.num_nodes() == 4);
    REQUIRE(prob.capacity() == 7.0);
}

// =====================================================================
// SEC separator: target_in_S = true case (path target inside cut set)
// =====================================================================

TEST_CASE("SEC path: rhs_coeff=1 when target in S (target_in_S=true)", "[sec][path]") {
    // source=0, target=3. Build a fractional solution where node 1 is
    // visited but has insufficient flow, and the min-cut for node 1
    // places {1,3} on the non-root side (S contains target=3).
    // With target in S, rhs_coeff should be 1 (not 2), so a flow of
    // 1.0 into S is enough and should NOT be violated.
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;  // source
    y_values[1] = 1.0;  // intermediate
    y_values[3] = 1.0;  // target
    // Node 2 not visited

    // Path: 0 -> 1 -> 3. Edges: {0,1} and {1,3} active.
    // Min-cut for node 1: S={1,3} has cut δ(S) = {edge(0,1)} with flow=1.
    // Since target=3 is in S, rhs_coeff = 1, violation = 1*1 - 1 = 0. No cut.
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 1 && v == 3)) {
            x_values[e] = 1.0;
        }
    }

    auto support = build_support(prob, x_values);

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // If rhs_coeff were incorrectly 2.0, we'd get violation 2*1-1=1 > 0.
    // Correct rhs_coeff=1.0 gives violation 0.0 → no cuts.
    REQUIRE(cuts.empty());
}

TEST_CASE("SEC path: target_in_S=false gives rhs_coeff=2", "[sec][path]") {
    // source=0, target=3. Build a solution where S does NOT contain target=3.
    // Node 2 is visited but disconnected; the min-cut for node 2 should give
    // S = {2} (not containing target=3), so rhs_coeff = 2.
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;  // source
    y_values[2] = 0.8;  // intermediate, not connected
    y_values[3] = 1.0;  // target

    // Only connect source to target directly: edge {0,3}
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 3) x_values[e] = 1.0;
    }

    auto support = build_support(prob, x_values);

    rcspp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    rcspp::sep::SECSeparator sec;
    auto cuts = sec.separate(ctx);

    // Node 2 is disconnected from source with y=0.8. S={2} doesn't contain
    // target=3, so rhs_coeff=2.0. Flow=0, violation = 2*0.8 - 0 = 1.6.
    REQUIRE(!cuts.empty());
    // Verify the violation magnitude is consistent with rhs_coeff=2
    REQUIRE(cuts[0].violation > 1.0);
}

// =====================================================================
// Preprocessing: demand reachability
// =====================================================================

TEST_CASE("Demand reachability: tour round-trip check", "[preprocess]") {
    // 3 nodes, depot=0. Demands: 0, 3, 5. Capacity=7.
    // Round-trip to node 1: 2*3 - 3 = 3 <= 7 → reachable
    // Round-trip to node 2: 2*5 - 5 = 5 <= 7 → reachable
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 3, 5};
    prob.build(3, edges, costs, profits, demands, 7.0, 0, 0);

    auto reachable = rcspp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);
    REQUIRE(reachable[1]);
    REQUIRE(reachable[2]);
}

TEST_CASE("Demand reachability: tour eliminates unreachable node", "[preprocess]") {
    // Capacity=4. Node 2 has demand=5.
    // Round-trip to node 2: 2*5 - 5 = 5 > 4 → unreachable
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 2, 5};
    prob.build(3, edges, costs, profits, demands, 4.0, 0, 0);

    auto reachable = rcspp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);
    REQUIRE(reachable[1]);  // 2*2 - 2 = 2 <= 4
    REQUIRE_FALSE(reachable[2]);  // 2*5 - 5 = 5 > 4
}

TEST_CASE("Demand reachability: path uses bidirectional check", "[preprocess][path]") {
    // 4 nodes, source=0, target=3. Demands: 0, 3, 4, 0. Capacity=6.
    // Path check: dist_s[v] + dist_t[v] - demand(v) <= Q
    //
    // With uniform edge demands = demand(endpoint), Dijkstra from source:
    // dist_s[0]=0, dist_s[1]=3, dist_s[2]=4, dist_s[3]=0 (demand=0)
    // Dijkstra from target (node 3):
    // dist_t[3]=0, dist_t[2]=4, dist_t[1]=3, dist_t[0]=0
    //
    // Node 1: dist_s=3, dist_t=3, check: 3+3-3=3 <= 6 → reachable
    // Node 2: dist_s=4, dist_t=4, check: 4+4-4=4 <= 6 → reachable
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    std::vector<double> costs = {1, 1, 1, 1, 1, 1};
    std::vector<double> profits = {0, 1, 1, 0};
    std::vector<double> demands = {0, 3, 4, 0};
    prob.build(4, edges, costs, profits, demands, 6.0, 0, 3);

    auto reachable = rcspp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);  // source always reachable
    REQUIRE(reachable[3]);  // target always reachable
    REQUIRE(reachable[1]);  // 3 + 3 - 3 = 3 <= 6
    REQUIRE(reachable[2]);  // 4 + 4 - 4 = 4 <= 6
}

TEST_CASE("Demand reachability: path eliminates unreachable node", "[preprocess][path]") {
    // Same setup but capacity=3. Node 2 has demand=4.
    // Node 2: dist_s=4, dist_t=4, check: 4+4-4=4 > 3 → unreachable
    // Node 1: dist_s=3, dist_t=3, check: 3+3-3=3 <= 3 → reachable
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    std::vector<double> costs = {1, 1, 1, 1, 1, 1};
    std::vector<double> profits = {0, 1, 1, 0};
    std::vector<double> demands = {0, 3, 4, 0};
    prob.build(4, edges, costs, profits, demands, 3.0, 0, 3);

    auto reachable = rcspp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);
    REQUIRE(reachable[3]);
    REQUIRE(reachable[1]);
    REQUIRE_FALSE(reachable[2]);
}

// =====================================================================
// Preprocessing: edge elimination
// =====================================================================

TEST_CASE("Edge elimination: tour eliminates expensive edge", "[preprocess]") {
    // 3 nodes. Edges: {0,1} cost=1, {0,2} cost=1, {1,2} cost=100.
    // Profits: all 0. Capacity: large.
    // With UB = 4 (tight), the edge {1,2} with cost=100 should be eliminated
    // since any tour using it costs at least 100 + something > 4.
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto eliminated = rcspp::preprocess::edge_elimination(prob, 4.0);
    REQUIRE_FALSE(eliminated[0]);  // {0,1} cost=1 — cheap
    REQUIRE_FALSE(eliminated[1]);  // {0,2} cost=1 — cheap
    REQUIRE(eliminated[2]);         // {1,2} cost=100 — eliminated
}

TEST_CASE("Edge elimination: path uses bidirectional labeling", "[preprocess][path]") {
    // 4 nodes, source=0, target=3.
    // Edges: {0,1}=1, {0,2}=1, {0,3}=1, {1,2}=100, {1,3}=1, {2,3}=1
    // With tight UB, edge {1,2} (cost=100) should be eliminated.
    rcspp::Problem prob;
    std::vector<rcspp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 3);

    auto eliminated = rcspp::preprocess::edge_elimination(prob, 5.0);
    // Edge {1,2} cost=100 should be eliminated — any path using it costs >= 100
    REQUIRE(eliminated[3]);
    // Cheap edges should survive
    REQUIRE_FALSE(eliminated[0]);  // {0,1}
    REQUIRE_FALSE(eliminated[2]);  // {0,3}
}

// =====================================================================
// Warm-start heuristic
// =====================================================================

TEST_CASE("Warm-start: tour produces closed loop", "[heuristic]") {
    auto prob = make_small_problem();
    auto result = rcspp::heuristic::build_warm_start(prob, 50);

    // Should produce a valid solution
    REQUIRE(result.objective < std::numeric_limits<double>::max());
    REQUIRE(result.col_values.size() ==
            static_cast<size_t>(prob.num_edges() + prob.num_nodes()));

    // Source/depot y-variable should be fixed
    REQUIRE(result.col_values[prob.num_edges() + prob.source()] == 1.0);
}

TEST_CASE("Warm-start: path produces valid open path", "[heuristic][path]") {
    auto prob = make_small_path_problem();
    auto result = rcspp::heuristic::build_warm_start(prob, 50);

    REQUIRE(result.objective < std::numeric_limits<double>::max());

    int32_t m = prob.num_edges();
    // Both source and target y-variables should be 1
    REQUIRE(result.col_values[m + prob.source()] == 1.0);
    REQUIRE(result.col_values[m + prob.target()] == 1.0);

    // Count active edges and visited nodes
    int active_edges = 0;
    for (int32_t e = 0; e < m; ++e) {
        if (result.col_values[e] > 0.5) active_edges++;
    }
    int visited_nodes = 0;
    for (int32_t i = 0; i < prob.num_nodes(); ++i) {
        if (result.col_values[m + i] > 0.5) visited_nodes++;
    }
    // Path: #edges = #visited_nodes - 1 (open path, no closing edge)
    REQUIRE(active_edges == visited_nodes - 1);
}

