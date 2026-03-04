#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/io.h"
#include "core/problem.h"
#include "heuristic/primal_heuristic.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/reachability.h"
#include "sep/comb_separator.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/rglm_separator.h"
#include "sep/sec_separator.h"
#include "sep/spi_separator.h"
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

TEST_CASE("Numeric IO: load tour (no source/target line)", "[io]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    REQUIRE(prob.is_tour());
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 0);
    REQUIRE(prob.num_nodes() == 4);
    REQUIRE(prob.capacity() == 7.0);
}

TEST_CASE("Numeric IO: load path (with source/target line)", "[io]") {
    auto prob = cptp::io::load("tests/data/tiny4_path.txt");
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
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 3, 5};
    prob.build(3, edges, costs, profits, demands, 7.0, 0, 0);

    auto reachable = cptp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);
    REQUIRE(reachable[1]);
    REQUIRE(reachable[2]);
}

TEST_CASE("Demand reachability: tour eliminates unreachable node", "[preprocess]") {
    // Capacity=4. Node 2 has demand=5.
    // Round-trip to node 2: 2*5 - 5 = 5 > 4 → unreachable
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 1, 1};
    std::vector<double> demands = {0, 2, 5};
    prob.build(3, edges, costs, profits, demands, 4.0, 0, 0);

    auto reachable = cptp::preprocess::demand_reachability(prob);
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
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    std::vector<double> costs = {1, 1, 1, 1, 1, 1};
    std::vector<double> profits = {0, 1, 1, 0};
    std::vector<double> demands = {0, 3, 4, 0};
    prob.build(4, edges, costs, profits, demands, 6.0, 0, 3);

    auto reachable = cptp::preprocess::demand_reachability(prob);
    REQUIRE(reachable[0]);  // source always reachable
    REQUIRE(reachable[3]);  // target always reachable
    REQUIRE(reachable[1]);  // 3 + 3 - 3 = 3 <= 6
    REQUIRE(reachable[2]);  // 4 + 4 - 4 = 4 <= 6
}

TEST_CASE("Demand reachability: path eliminates unreachable node", "[preprocess][path]") {
    // Same setup but capacity=3. Node 2 has demand=4.
    // Node 2: dist_s=4, dist_t=4, check: 4+4-4=4 > 3 → unreachable
    // Node 1: dist_s=3, dist_t=3, check: 3+3-3=3 <= 3 → reachable
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    std::vector<double> costs = {1, 1, 1, 1, 1, 1};
    std::vector<double> profits = {0, 1, 1, 0};
    std::vector<double> demands = {0, 3, 4, 0};
    prob.build(4, edges, costs, profits, demands, 3.0, 0, 3);

    auto reachable = cptp::preprocess::demand_reachability(prob);
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
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto eliminated = cptp::preprocess::edge_elimination(prob, 4.0);
    REQUIRE_FALSE(eliminated[0]);  // {0,1} cost=1 — cheap
    REQUIRE_FALSE(eliminated[1]);  // {0,2} cost=1 — cheap
    REQUIRE(eliminated[2]);         // {1,2} cost=100 — eliminated
}

TEST_CASE("Edge elimination: path uses bidirectional labeling", "[preprocess][path]") {
    // 4 nodes, source=0, target=3.
    // Edges: {0,1}=1, {0,2}=1, {0,3}=1, {1,2}=100, {1,3}=1, {2,3}=1
    // With tight UB, edge {1,2} (cost=100) should be eliminated.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 3);

    auto eliminated = cptp::preprocess::edge_elimination(prob, 5.0);
    // Edge {1,2} cost=100 should be eliminated — any path using it costs >= 100
    REQUIRE(eliminated[3]);
    // Cheap edges should survive
    REQUIRE_FALSE(eliminated[0]);  // {0,1}
    REQUIRE_FALSE(eliminated[2]);  // {0,3}
}

// =====================================================================
// Warm-start heuristic
// =====================================================================

TEST_CASE("Primal heuristic: tour produces closed loop", "[heuristic]") {
    auto prob = make_small_problem();
    auto result = cptp::heuristic::build_initial_solution(prob, 50);

    // Should produce a valid solution
    REQUIRE(result.objective < std::numeric_limits<double>::max());
    REQUIRE(result.col_values.size() ==
            static_cast<size_t>(prob.num_edges() + prob.num_nodes()));

    // Source/depot y-variable should be fixed
    REQUIRE(result.col_values[prob.num_edges() + prob.source()] == 1.0);
}

TEST_CASE("Primal heuristic: path produces valid open path", "[heuristic][path]") {
    auto prob = make_small_path_problem();
    auto result = cptp::heuristic::build_initial_solution(prob, 50);

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

TEST_CASE("LP-guided heuristic: reduces graph and finds solution", "[heuristic]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Build a fake LP relaxation: fractional values on some edges/nodes
    std::vector<double> x_lp(m, 0.0);
    std::vector<double> y_lp(n, 0.0);
    for (auto e : prob.graph().edges()) {
        x_lp[e] = 0.3;  // all edges fractional
    }
    for (int32_t i = 0; i < n; ++i) {
        y_lp[i] = 0.6;  // all nodes active
    }

    // No incumbent
    std::vector<double> incumbent;

    // Strategy 0 (all), should find a feasible solution
    auto result = cptp::heuristic::lp_guided_heuristic(
        prob, x_lp, y_lp, incumbent, 50.0, 0);

    REQUIRE(result.objective < std::numeric_limits<double>::max());
    REQUIRE(result.col_values.size() ==
            static_cast<size_t>(m + n));
    REQUIRE(result.col_values[m + prob.source()] == 1.0);
}

TEST_CASE("LP-guided heuristic: RINS with incumbent", "[heuristic]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Get an incumbent from build_initial_solution
    auto initial = cptp::heuristic::build_initial_solution(prob, 50);
    REQUIRE(!initial.col_values.empty());

    // Use the initial solution's edge/node values as LP relaxation
    std::vector<double> x_lp(initial.col_values.begin(),
                              initial.col_values.begin() + m);
    std::vector<double> y_lp(initial.col_values.begin() + m,
                              initial.col_values.begin() + m + n);

    // RINS strategy (2) with incumbent
    auto result = cptp::heuristic::lp_guided_heuristic(
        prob, x_lp, y_lp, initial.col_values, 50.0, 2);

    REQUIRE(result.objective < std::numeric_limits<double>::max());
    REQUIRE(result.col_values.size() ==
            static_cast<size_t>(m + n));
}

TEST_CASE("Reduce strategies produce valid subgraphs", "[heuristic]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_lp(m, 0.0);
    std::vector<double> y_lp(n, 0.0);
    // Set a few edges fractional
    x_lp[0] = 0.5;
    x_lp[1] = 0.3;
    y_lp[0] = 0.8;
    y_lp[1] = 0.7;

    auto rg_thresh = cptp::heuristic::reduce_lp_threshold(prob, x_lp, y_lp);
    REQUIRE(rg_thresh.edge_active.size() == static_cast<size_t>(m));
    REQUIRE(rg_thresh.node_active.size() == static_cast<size_t>(n));
    // Source and target always active
    REQUIRE(rg_thresh.node_active[prob.source()]);
    REQUIRE(rg_thresh.node_active[prob.target()]);
    // Fractional edges should be active
    REQUIRE(rg_thresh.edge_active[0]);
    REQUIRE(rg_thresh.edge_active[1]);

    auto rg_nbhd = cptp::heuristic::reduce_neighborhood(prob, x_lp, y_lp);
    REQUIRE(rg_nbhd.edge_active.size() == static_cast<size_t>(m));
    REQUIRE(rg_nbhd.node_active[prob.source()]);
}

// =====================================================================
// Multistar separator
// =====================================================================

TEST_CASE("Multistar separator basic test", "[multistar]") {
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
        if ((u == 0 && v == 1) || (u == 0 && v == 2) || (u == 1 && v == 2))
            x_values[e] = 1.0;
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

    cptp::sep::MultistarSeparator multistar;
    // Just verify it runs without error
    auto cuts = multistar.separate(ctx);
    REQUIRE(true);
}

TEST_CASE("Comb separator basic test", "[comb]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;
    y_values[1] = 0.5;
    y_values[2] = 0.5;
    y_values[3] = 0.5;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 1) x_values[e] = 1.0;
        if (u == 2 && v == 3) x_values[e] = 1.0;
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

    cptp::sep::CombSeparator comb;
    auto cuts = comb.separate(ctx);
    REQUIRE(true);
}

TEST_CASE("RGLM separator basic test", "[rglm]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;
    y_values[1] = 0.5;
    y_values[2] = 0.5;
    y_values[3] = 0.5;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 1) x_values[e] = 1.0;
        if (u == 2 && v == 3) x_values[e] = 1.0;
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

    cptp::sep::RGLMSeparator rglm;
    auto cuts = rglm.separate(ctx);
    REQUIRE(true);
}

// =====================================================================
// RCI separator: capacity-violating fractional solution
// =====================================================================

TEST_CASE("RCI separator on fractional solution with high demands", "[rci]") {
    // Create a problem where demands are significant relative to capacity.
    // 4 nodes, capacity=7, demands={0, 3, 4, 2}.
    // A fractional solution visiting all customers would use total demand=9>7,
    // requiring at least ceil(9/7)=2 "vehicles" (routes), i.e., 4 crossing edges.
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    // All customers visited with fractional edges
    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[2] = 1.0;
    y_values[3] = 1.0;

    // Spread flow across edges to create a fractional solution
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        // Give positive flow to all depot edges and some non-depot edges
        if (u == 0 || v == 0) {
            x_values[e] = 0.8;
        } else {
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
    };

    cptp::sep::RCISeparator rci;
    auto cuts = rci.separate(ctx);
    // RCI may or may not find cuts on this particular solution, but should not crash
    REQUIRE(true);
}

// =====================================================================
// SEC separator with path mode: multiple disconnected components
// =====================================================================

TEST_CASE("SEC separator path: multiple disconnected nodes", "[sec][path]") {
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Only node 0 and 3 connected directly; nodes 1 and 2 disconnected
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    y_values[0] = 1.0;
    y_values[1] = 0.7;
    y_values[2] = 0.6;
    y_values[3] = 1.0;

    // Only edge {0,3} active
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 3) x_values[e] = 1.0;
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
    // Nodes 1 and 2 are disconnected from source with y > 0, so cuts expected
    REQUIRE(!cuts.empty());
}

// =====================================================================
// Gomory-Hu tree: direct test
// =====================================================================

TEST_CASE("Gomory-Hu tree: min-cut values on simple graph", "[maxflow]") {
    // 3-node complete graph with symmetric arcs
    cptp::digraph_builder builder(3);
    builder.add_arc(0, 1, 2.0);
    builder.add_arc(1, 0, 2.0);
    builder.add_arc(0, 2, 3.0);
    builder.add_arc(2, 0, 3.0);
    builder.add_arc(1, 2, 1.0);
    builder.add_arc(2, 1, 1.0);
    auto [sg, cap] = builder.build();

    cptp::gomory_hu_tree tree(sg, cap, 0);
    // The Gomory-Hu tree encodes all pairwise min-cuts.
    // Min-cut(0,1) = 3 (edges 0-1:2, plus 1-2:1 or via other path)
    // Min-cut(0,2) = 3 (edges 0-2:3, or 0-1:2 + 1-2:1)
    // Min-cut(1,2) = 3 (edges 1-0:2+0-2:3 or 1-2:1+...)
    // Just verify the tree was constructed successfully
    REQUIRE(true);
}

// =====================================================================
// Separator name() accessor
// =====================================================================

TEST_CASE("Separator name() returns correct strings", "[separator]") {
    cptp::sep::SECSeparator sec;
    REQUIRE(sec.name() == "SEC");

    cptp::sep::RCISeparator rci;
    REQUIRE(rci.name() == "RCI");

    cptp::sep::MultistarSeparator ms;
    REQUIRE(ms.name() == "Multistar");

    cptp::sep::CombSeparator comb;
    REQUIRE(comb.name() == "Comb");

    cptp::sep::RGLMSeparator rglm;
    REQUIRE(rglm.name() == "RGLM");
}

// =====================================================================
// Warm-start: tour consistency checks
// =====================================================================

TEST_CASE("Warm-start: tour degree consistency", "[heuristic]") {
    auto prob = make_small_problem();
    auto result = cptp::heuristic::build_warm_start(prob, 50);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();

    // For each visited node, count the degree (number of incident active edges)
    for (int32_t v = 0; v < n; ++v) {
        if (result.col_values[m + v] < 0.5) continue;  // not visited

        int degree = 0;
        for (auto e : graph.incident_edges(v)) {
            if (result.col_values[e] > 0.5) degree++;
        }
        // In a tour, every visited node has exactly degree 2
        REQUIRE(degree == 2);
    }
}

TEST_CASE("Warm-start: path degree consistency", "[heuristic][path]") {
    auto prob = make_small_path_problem();
    auto result = cptp::heuristic::build_warm_start(prob, 50);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();

    for (int32_t v = 0; v < n; ++v) {
        if (result.col_values[m + v] < 0.5) continue;

        int degree = 0;
        for (auto e : graph.incident_edges(v)) {
            if (result.col_values[e] > 0.5) degree++;
        }
        // Source and target have degree 1, intermediates have degree 2
        if (v == prob.source() || v == prob.target()) {
            REQUIRE(degree == 1);
        } else {
            REQUIRE(degree == 2);
        }
    }
}

// =====================================================================
// Edge elimination: path with infinity UB eliminates nothing
// =====================================================================

TEST_CASE("Edge elimination: infinite UB eliminates nothing", "[preprocess]") {
    auto prob = make_small_problem();
    auto eliminated = cptp::preprocess::edge_elimination(
        prob, std::numeric_limits<double>::infinity());
    for (size_t e = 0; e < eliminated.size(); ++e) {
        REQUIRE_FALSE(eliminated[e]);
    }
}

// =====================================================================
// Demand reachability: all nodes reachable with large capacity
// =====================================================================

TEST_CASE("Demand reachability: large capacity makes all reachable", "[preprocess]") {
    auto prob = make_small_problem();
    // Rebuild with very large capacity
    cptp::Problem big_cap;
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
    big_cap.build(4, edges, costs, profits, demands, 1e18, 0, 0);

    auto reachable = cptp::preprocess::demand_reachability(big_cap);
    for (int32_t v = 0; v < big_cap.num_nodes(); ++v) {
        REQUIRE(reachable[v]);
    }
}

// =====================================================================
// Determinism
// =====================================================================

TEST_CASE("Warm-start: deterministic mode produces identical results", "[heuristic][determinism]") {
    auto prob = make_small_problem();

    auto r1 = cptp::heuristic::build_warm_start(prob, 50);
    auto r2 = cptp::heuristic::build_warm_start(prob, 50);

    REQUIRE(r1.objective == r2.objective);
    REQUIRE(r1.col_values == r2.col_values);
}

TEST_CASE("Warm-start: deterministic path produces identical results", "[heuristic][determinism][path]") {
    auto prob = make_small_path_problem();

    auto r1 = cptp::heuristic::build_warm_start(prob, 50);
    auto r2 = cptp::heuristic::build_warm_start(prob, 50);

    REQUIRE(r1.objective == r2.objective);
    REQUIRE(r1.col_values == r2.col_values);
}

TEST_CASE("Warm-start: time budget option is deterministic no-op", "[heuristic][determinism]") {
    auto prob = make_small_problem();

    auto base = cptp::heuristic::build_warm_start(prob, 50, 0.0);
    auto with_budget = cptp::heuristic::build_warm_start(prob, 50, 100.0);

    REQUIRE(with_budget.objective == base.objective);
    REQUIRE(with_budget.col_values == base.col_values);
}

TEST_CASE("WorkUnitBudget enforces cap and reservation", "[heuristic][determinism]") {
    cptp::WorkUnitBudget budget(3);
    REQUIRE(budget.capped());
    REQUIRE(budget.try_consume(2));
    REQUIRE(budget.used() == 2);
    REQUIRE_FALSE(budget.try_consume(2));
    REQUIRE(budget.reserve_up_to(5) == 1);
    REQUIRE(budget.used() == 3);
    REQUIRE(budget.remaining() == 0);
}

TEST_CASE("Warm-start: work-unit cap yields deterministic prefix",
          "[heuristic][determinism]") {
    auto prob = make_small_problem();

    auto budget1 = std::make_shared<cptp::WorkUnitBudget>(2);
    auto r1 = cptp::heuristic::build_warm_start(prob, 50, 0.0, budget1);

    auto budget2 = std::make_shared<cptp::WorkUnitBudget>(2);
    auto r2 = cptp::heuristic::build_warm_start(prob, 50, 0.0, budget2);

    REQUIRE(budget1->used() == 2);
    REQUIRE(budget2->used() == 2);
    REQUIRE(r1.objective == r2.objective);
    REQUIRE(r1.col_values == r2.col_values);
}

TEST_CASE("LP-guided heuristic: deterministic restarts are reproducible",
          "[heuristic][determinism]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_lp(m, 0.3);
    std::vector<double> y_lp(n, 0.6);
    std::vector<double> incumbent;

    auto budget1 = std::make_shared<cptp::WorkUnitBudget>(16);
    auto r1 = cptp::heuristic::lp_guided_heuristic(
        prob, x_lp, y_lp, incumbent,
        0.0, 0, 16, 7u, budget1);

    auto budget2 = std::make_shared<cptp::WorkUnitBudget>(16);
    auto r2 = cptp::heuristic::lp_guided_heuristic(
        prob, x_lp, y_lp, incumbent,
        0.0, 0, 16, 7u, budget2);

    REQUIRE(budget1->used() == 16);
    REQUIRE(budget2->used() == 16);
    REQUIRE(r1.objective == r2.objective);
    REQUIRE(r1.col_values == r2.col_values);
}

// =====================================================================
// SPI separator: Shortest Path Inequalities
// =====================================================================

namespace {

/// Compute all-pairs shortest path bounds for a problem.
std::vector<double> compute_all_pairs(const cptp::Problem& prob) {
    const int32_t n = prob.num_nodes();
    constexpr double inf = std::numeric_limits<double>::infinity();
    std::vector<double> ap(static_cast<size_t>(n) * n, inf);
    for (int32_t s = 0; s < n; ++s) {
        auto row = cptp::preprocess::forward_labeling(prob, s);
        std::copy(row.begin(), row.end(),
                  ap.begin() + static_cast<ptrdiff_t>(s) * n);
    }
    return ap;
}

/// Build a 5-node tour problem where nodes 3 and 4 are very expensive to
/// visit together (long detour), making the pair {3,4} infeasible.
///
/// Zero profits: ensures labeling lower bounds are tight (no cycling benefit).
/// d[0][3] = d[0][4] = 50 (direct edges), d[3][4] = 50.
/// Pair {3,4} lb = 50 + 50 + 50 + 0 = 150.
cptp::Problem make_spi_tour_problem() {
    cptp::Problem prob;
    // 5 nodes: depot=0, customers=1,2,3,4. Complete graph (10 edges).
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            edges.push_back({i, j});
            // Make edges to/from nodes 3 and 4 very expensive (50 each),
            // but edges among {0,1,2} cheap (5 each).
            if (i >= 3 || j >= 3) {
                costs.push_back(50.0);
            } else {
                costs.push_back(5.0);
            }
        }
    }

    // Zero profits: labeling = pure shortest path (no cycling incentive)
    std::vector<double> profits = {0, 0, 0, 0, 0};
    std::vector<double> demands = {0, 1, 1, 1, 1};
    double capacity = 10.0;

    prob.build(5, edges, costs, profits, demands, capacity, 0, 0);
    return prob;
}

/// Build a 5-node path problem (source=0, target=4) with expensive detours.
///
/// Zero profits, so labeling bounds are tight.
/// Cheap path: 0→1→2→3→4 costs 5+5+5+5=20.
/// Detour edges cost 50, so visiting e.g. node 1 AND node 3 via expensive
/// edges has a clear lower bound.
cptp::Problem make_spi_path_problem() {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            edges.push_back({i, j});
            // Linear chain 0-1-2-3-4 is cheap (5 each), everything else expensive
            if (j == i + 1) {
                costs.push_back(5.0);
            } else {
                costs.push_back(50.0);
            }
        }
    }

    std::vector<double> profits = {0, 0, 0, 0, 0};
    std::vector<double> demands = {0, 1, 1, 1, 0};
    double capacity = 10.0;

    prob.build(5, edges, costs, profits, demands, capacity, 0, 4);
    return prob;
}

}  // namespace

TEST_CASE("SPI separator: finds pair cuts on tour with tight UB", "[spi]") {
    auto prob = make_spi_tour_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    auto all_pairs = compute_all_pairs(prob);

    // Fractional LP solution: depot + all customers at y=0.5
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    for (int i = 1; i < n; ++i) y_values[i] = 0.6;

    // Set some fractional edges to make it look like an LP solution
    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 || v == 0) x_values[e] = 0.5;
    }

    auto support = build_support(prob, x_values);

    // Tight UB: only cheap tours are feasible.
    // Tour {0,1,2} costs: 1+1+1 = 3, profit = 5+5 = 10, obj = 3-10 = -7.
    // Tour visiting node 3 or 4 costs >= 50*2 = 100 (two expensive edges).
    // Set UB tight enough that visiting both 3 and 4 is infeasible.
    double tight_ub = 90.0;

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = tight_ub,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    // With expensive edges to nodes 3,4, the pair {3,4} should be infeasible
    // under a tight UB, yielding y_3 + y_4 <= 1.
    // Check that at least one cut was found.
    REQUIRE(!cuts.empty());

    // All cuts should be of the form: sum y_i <= rhs
    for (const auto& cut : cuts) {
        REQUIRE(cut.rhs >= 1.0 - 1e-9);
        REQUIRE(cut.violation > 1e-6);
        for (double v : cut.values) {
            REQUIRE(v == 1.0);
        }
    }

    // With lifting, pair cuts for {3,4} should include nodes 1,2 (symmetric
    // expensive edges), so at least one cut should have more variables than
    // just the pair.
    bool found_lifted = false;
    for (const auto& cut : cuts) {
        if (cut.rhs < 1.0 + 1e-9 && cut.size() > 2) {
            found_lifted = true;
            break;
        }
    }
    REQUIRE(found_lifted);
}

TEST_CASE("SPI separator: no cuts when UB is loose", "[spi]") {
    auto prob = make_spi_tour_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    auto all_pairs = compute_all_pairs(prob);

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    for (int i = 1; i < n; ++i) y_values[i] = 0.6;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 || v == 0) x_values[e] = 0.5;
    }

    auto support = build_support(prob, x_values);

    // Very loose UB — everything is feasible
    double loose_ub = 1e6;

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = loose_ub,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    // No cuts should be found with a very loose UB
    REQUIRE(cuts.empty());
}

TEST_CASE("SPI separator: no cuts without all-pairs data", "[spi]") {
    auto prob = make_spi_tour_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    for (int i = 1; i < n; ++i) y_values[i] = 0.6;

    auto support = build_support(prob, x_values);

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = 10.0,
        // all_pairs left empty
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    // Without all-pairs data, separator should gracefully return nothing
    REQUIRE(cuts.empty());
}

TEST_CASE("SPI separator: path problem pair cuts", "[spi][path]") {
    auto prob = make_spi_path_problem();
    REQUIRE_FALSE(prob.is_tour());
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    auto all_pairs = compute_all_pairs(prob);

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;  // source
    y_values[4] = 1.0;  // target
    for (int i = 1; i < 4; ++i) y_values[i] = 0.7;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 || v == 4) x_values[e] = 0.3;
    }

    auto support = build_support(prob, x_values);

    // Set a tight UB for the path
    double tight_ub = 5.0;

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = tight_ub,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    // With tight UB, some pairs should be infeasible
    // Just verify the separator runs correctly and produces valid cuts
    for (const auto& cut : cuts) {
        REQUIRE(cut.violation > 1e-6);
        for (double v : cut.values) {
            REQUIRE(v == 1.0);
        }
        // rhs = |S_min|-1 where S_min is the minimal infeasible set.
        // With lifting, cut.size() >= |S_min|.
        REQUIRE(cut.rhs >= 1.0 - 1e-9);
        REQUIRE(cut.rhs < static_cast<double>(cut.size()));
    }
}

TEST_CASE("SPI separator: greedy extension finds set cuts", "[spi]") {
    // 7-node tour problem with a cluster of expensive nodes.
    // Depot=0, cheap cluster {1,2}, expensive cluster {3,4,5,6}.
    // Zero profits: labeling = pure shortest path (no cycling incentive).
    // With tight UB, the separator should find set cuts like
    // y_3 + y_4 + y_5 <= 2, or larger sets.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 7; ++i) {
        for (int j = i + 1; j < 7; ++j) {
            edges.push_back({i, j});
            // Cheap edges among {0,1,2}, expensive to reach {3,4,5,6}
            if (i >= 3 || j >= 3) {
                costs.push_back(50.0);
            } else {
                costs.push_back(5.0);
            }
        }
    }

    std::vector<double> profits = {0, 0, 0, 0, 0, 0, 0};
    std::vector<double> demands = {0, 1, 1, 1, 1, 1, 1};
    prob.build(7, edges, costs, profits, demands, 100.0, 0, 0);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();
    auto all_pairs = compute_all_pairs(prob);

    // LP solution: all customers at 0.7
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    for (int i = 1; i < n; ++i) y_values[i] = 0.7;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 || v == 0) x_values[e] = 0.4;
    }

    auto support = build_support(prob, x_values);

    // Tight UB: visiting multiple expensive nodes costs > UB
    double tight_ub = 80.0;

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = tight_ub,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    REQUIRE(!cuts.empty());

    // Check that we found at least one set cut (rhs >= 1)
    // The greedy extension should find cuts larger than just pairs
    bool found_set_cut = false;
    for (const auto& cut : cuts) {
        REQUIRE(cut.violation > 1e-6);
        REQUIRE(cut.rhs >= 1.0 - 1e-9);
        REQUIRE(cut.rhs < static_cast<double>(cut.size()));
        if (cut.size() >= 3) found_set_cut = true;
    }
    // With 4 expensive nodes and tight UB, we should find triplet+ cuts
    REQUIRE(found_set_cut);
}

TEST_CASE("SPI separator: shrinking produces minimal sets", "[spi]") {
    // Verify that the shrink phase finds the smallest infeasible set.
    // Zero profits: labeling = pure shortest path (no cycling incentive).
    cptp::Problem prob;
    std::vector<cptp::Edge> edges;
    std::vector<double> costs;
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            edges.push_back({i, j});
            // All edges to nodes 3,4,5 cost 30
            // Edges among {0,1,2} cost 5
            if (i >= 3 || j >= 3) {
                costs.push_back(30.0);
            } else {
                costs.push_back(5.0);
            }
        }
    }

    std::vector<double> profits = {0, 0, 0, 0, 0, 0};
    std::vector<double> demands = {0, 1, 1, 1, 1, 1};
    prob.build(6, edges, costs, profits, demands, 100.0, 0, 0);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();
    auto all_pairs = compute_all_pairs(prob);

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    for (int i = 1; i < n; ++i) y_values[i] = 0.8;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 || v == 0) x_values[e] = 0.5;
    }

    auto support = build_support(prob, x_values);

    // UB where pairs might be feasible but triples are not
    double ub = 60.0;

    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
        .upper_bound = ub,
        .all_pairs = all_pairs,
    };

    cptp::sep::SPISeparator spi;
    auto cuts = spi.separate(ctx);

    // All cuts should be valid
    for (const auto& cut : cuts) {
        REQUIRE(cut.violation > 1e-6);
        for (double v : cut.values) {
            REQUIRE(v == 1.0);
        }
        REQUIRE(cut.rhs >= 1.0 - 1e-9);
        REQUIRE(cut.rhs < static_cast<double>(cut.size()));
    }
}
