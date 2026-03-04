#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <memory>

#include "core/problem.h"
#include "preprocess/bound_propagator.h"
#include "preprocess/edge_elimination.h"
#include "sep/sec_separator.h"
#include "sep/separation_oracle.h"

namespace {

/// Build a small 4-node tour problem: depot=0, customers=1,2,3.
cptp::Problem make_small_problem() {
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
    prob.build(4, edges, costs, profits, demands, 7.0, 0, 0);
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
    prob.build(4, edges, costs, profits, demands, 7.0, 0, 3);
    return prob;
}

}  // namespace

// =====================================================================
// SeparationOracle tests
// =====================================================================

TEST_CASE("SeparationOracle: finds cuts on violated solution", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();

    // Fractional solution with disconnected subtour
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

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    REQUIRE(!cuts.empty());

    // Cuts should be sorted by violation (descending)
    for (size_t i = 1; i < cuts.size(); ++i) {
        REQUIRE(cuts[i - 1].violation >= cuts[i].violation);
    }
}

TEST_CASE("SeparationOracle: no cuts on feasible solution", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();

    // Feasible tour: 0 - 1 - 2 - 0
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[2] = 1.0;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 0 && v == 2) || (u == 1 && v == 2))
            x_values[e] = 1.0;
    }

    auto cuts = oracle.separate(x_values, y_values, 0, m, 1e-6);
    REQUIRE(cuts.empty());
}

TEST_CASE("SeparationOracle: is_feasible works correctly", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);

    // Feasible tour: 0 - 1 - 2 - 0
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[2] = 1.0;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 0 && v == 2) || (u == 1 && v == 2))
            x_values[e] = 1.0;
    }

    REQUIRE(oracle.is_feasible(x_values, y_values, 0, m));

    // Now make it infeasible: disconnect node 2
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 1 && v == 2) x_values[e] = 0.0;
        if (u == 0 && v == 2) x_values[e] = 0.0;
        if (u == 2 && v == 3) x_values[e] = 1.0;
    }
    y_values[3] = 1.0;

    REQUIRE_FALSE(oracle.is_feasible(x_values, y_values, 0, m));
}

TEST_CASE("SeparationOracle: max_cuts_per_separator limits output", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();
    oracle.set_max_cuts_per_separator(1);

    // Fractional solution with violations
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

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    // With 4 separators at max 1 each, should get at most 4 cuts
    REQUIRE(cuts.size() <= 4);
}

TEST_CASE("SeparationOracle: path mode works", "[oracle][path]") {
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();
    REQUIRE(oracle.separators().size() == 4);

    // Disconnected solution — node 2 not connected to source
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[3] = 1.0;
    y_values[1] = 0.5;
    y_values[2] = 0.5;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 3) || (u == 1 && v == 2))
            x_values[e] = 1.0;
    }

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    REQUIRE(!cuts.empty());
}

// =====================================================================
// BoundPropagator tests
// =====================================================================

TEST_CASE("BoundPropagator: sweep fixes expensive edges", "[propagator]") {
    // 3 nodes, tour. Edges: {0,1}=1, {0,2}=1, {1,2}=100.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);
    auto bwd = fwd;  // tour: same
    double correction = prob.profit(0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, bwd, correction);

    // col_upper for edges (all unfixed)
    std::vector<double> col_upper = {1.0, 1.0, 1.0};
    double upper_bound = 4.0;

    auto fixings = prop.sweep(upper_bound, col_upper);

    // Edge {1,2} with cost=100 should be fixed to 0
    bool fixes_edge_2 = false;
    for (int32_t e : fixings) {
        if (e == 2) fixes_edge_2 = true;
    }
    REQUIRE(fixes_edge_2);

    // Cheap edges should NOT be fixed
    bool fixes_edge_0 = false;
    bool fixes_edge_1 = false;
    for (int32_t e : fixings) {
        if (e == 0) fixes_edge_0 = true;
        if (e == 1) fixes_edge_1 = true;
    }
    REQUIRE_FALSE(fixes_edge_0);
    REQUIRE_FALSE(fixes_edge_1);
}

TEST_CASE("BoundPropagator: sweep_nodes fixes isolated nodes", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    // Simulate: edge 0 and edge 2 are fixed to 0, edge 1 is still open.
    // Node 1 is incident to edges 0 ({0,1}) and 2 ({1,2}).
    // With both fixed, node 1 should be fixable.
    int32_t m = prob.num_edges();
    std::vector<double> col_upper(m + prob.num_nodes(), 1.0);
    col_upper[0] = 0.0;  // edge {0,1} fixed
    col_upper[2] = 0.0;  // edge {1,2} fixed

    auto node_fixings = prop.sweep_nodes(col_upper, m);

    // Node 1 (col index m+1) should be in fixings
    bool fixes_node_1 = false;
    for (int32_t idx : node_fixings) {
        if (idx == m + 1) fixes_node_1 = true;
    }
    REQUIRE(fixes_node_1);
}

TEST_CASE("BoundPropagator: propagate_fixed_edge (neighbor-only)", "[propagator]") {
    // 4 nodes, tour. Known expensive path through certain edges.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    // All edges unfixed
    int32_t m = prob.num_edges();
    std::vector<double> col_upper(m, 1.0);

    // Fix edge 0 ({0,1}) to 1. With tight UB, edge 3 ({1,2} cost=100)
    // should be fixable because any path using both {0,1} and {1,2} is expensive.
    double upper_bound = 5.0;
    auto fixings = prop.propagate_fixed_edge(0, upper_bound, col_upper);

    // Edge {1,2} (index 3) costs 100 and should be eliminated
    bool fixes_expensive = false;
    for (int32_t e : fixings) {
        if (e == 3) fixes_expensive = true;
    }
    REQUIRE(fixes_expensive);
}

TEST_CASE("BoundPropagator: propagate_fixed_edge (all-pairs)", "[propagator]") {
    // 4 nodes, tour. Same setup but with all-pairs bounds.
    // The all-pairs path scans ALL edges (not just neighbors), so it should
    // also eliminate the expensive edge {1,2} even though it's not adjacent
    // to the fixed edge {0,3}.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    // Make {1,2} very expensive — any tour through it costs >= 100
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 0);

    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    // Build all-pairs bounds: run forward labeling from each node
    std::vector<double> all_pairs(n * n, std::numeric_limits<double>::max());
    for (int32_t s = 0; s < n; ++s) {
        auto bounds = cptp::preprocess::forward_labeling(prob, s);
        for (int32_t v = 0; v < n; ++v) {
            all_pairs[s * n + v] = bounds[v];
        }
    }

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);
    prop.set_all_pairs_bounds(all_pairs);
    REQUIRE(prop.has_all_pairs_bounds());

    std::vector<double> col_upper(m, 1.0);

    // Fix edge 2 ({0,3}) to 1. With tight UB=5, any tour through {0,3}
    // AND {1,2} costs at least 100+... > 5, so {1,2} should be eliminated.
    // This tests the all-pairs code path (edge {1,2} is NOT adjacent to {0,3}).
    double upper_bound = 5.0;
    auto fixings = prop.propagate_fixed_edge(2, upper_bound, col_upper);

    bool fixes_expensive = false;
    for (int32_t e : fixings) {
        if (e == 3) fixes_expensive = true;  // edge {1,2}
    }
    REQUIRE(fixes_expensive);
}

// =====================================================================
// SeparationOracle: accumulation behavior
// =====================================================================

TEST_CASE("SeparationOracle: add_default_separators is cumulative", "[oracle]") {
    auto prob = make_small_problem();

    cptp::sep::SeparationOracle oracle(prob);
    REQUIRE(oracle.separators().size() == 0);

    oracle.add_default_separators();
    REQUIRE(oracle.separators().size() == 4);

    // Calling again adds another 4
    oracle.add_default_separators();
    REQUIRE(oracle.separators().size() == 8);
}

TEST_CASE("SeparationOracle: individual separator addition", "[oracle]") {
    auto prob = make_small_problem();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_separator(std::make_unique<cptp::sep::SECSeparator>());
    REQUIRE(oracle.separators().size() == 1);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    // Fractional solution with disconnected subtour — SEC alone should find cuts
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

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    REQUIRE(!cuts.empty());
}

// =====================================================================
// SeparationOracle: edge cases and additional coverage
// =====================================================================

TEST_CASE("SeparationOracle: no separators returns no cuts", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    // Do NOT add any separators

    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[1] = 0.5;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 0 && v == 1) x_values[e] = 1.0;
        if (u == 2 && v == 3) x_values[e] = 1.0;
    }

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    REQUIRE(cuts.empty());
}

TEST_CASE("SeparationOracle: cut struct has valid fields", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();

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

    auto cuts = oracle.separate(x_values, y_values, 0, m);
    REQUIRE(!cuts.empty());

    for (const auto& cut : cuts) {
        // indices and values must have equal length
        REQUIRE(cut.indices.size() == cut.values.size());
        REQUIRE(cut.size() == static_cast<int32_t>(cut.indices.size()));
        // non-empty cut
        REQUIRE(cut.size() > 0);
        // violation should be positive
        REQUIRE(cut.violation > 0.0);
        // all indices should be in valid range [0, m+n)
        for (int32_t idx : cut.indices) {
            REQUIRE(idx >= 0);
            REQUIRE(idx < m + n);
        }
    }
}

TEST_CASE("SeparationOracle: is_feasible works in path mode", "[oracle][path]") {
    auto prob = make_small_path_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);

    // Feasible path: 0 -> 1 -> 3
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);
    y_values[0] = 1.0;
    y_values[1] = 1.0;
    y_values[3] = 1.0;

    const auto& graph = prob.graph();
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if ((u == 0 && v == 1) || (u == 1 && v == 3))
            x_values[e] = 1.0;
    }

    REQUIRE(oracle.is_feasible(x_values, y_values, 0, m));

    // Make infeasible: disconnect node 2 from source, activate it
    for (auto e : graph.edges()) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        if (u == 1 && v == 3) x_values[e] = 0.0;
        if (u == 1 && v == 2) x_values[e] = 1.0;
    }
    y_values[1] = 0.0;
    y_values[2] = 1.0;

    REQUIRE_FALSE(oracle.is_feasible(x_values, y_values, 0, m));
}

TEST_CASE("SeparationOracle: separate with non-zero x_offset", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();

    // Fractional solution with disconnected subtour
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

    // Use non-zero offsets
    int32_t x_off = 10;
    int32_t y_off = x_off + m;
    auto cuts = oracle.separate(x_values, y_values, x_off, y_off);
    REQUIRE(!cuts.empty());

    // All cut indices should be offset accordingly
    for (const auto& cut : cuts) {
        for (int32_t idx : cut.indices) {
            REQUIRE(idx >= x_off);
        }
    }
}

TEST_CASE("SeparationOracle: max_cuts_per_separator 0 means unlimited", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();
    oracle.set_max_cuts_per_separator(0);
    REQUIRE(oracle.max_cuts_per_separator() == 0);

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

    // Should not crash, and may return more cuts than with limit=1
    auto cuts_unlimited = oracle.separate(x_values, y_values, 0, m);

    cptp::sep::SeparationOracle oracle2(prob);
    oracle2.add_default_separators();
    oracle2.set_max_cuts_per_separator(1);
    auto cuts_limited = oracle2.separate(x_values, y_values, 0, m);

    REQUIRE(cuts_unlimited.size() >= cuts_limited.size());
}

TEST_CASE("SeparationOracle: all-zero solution produces no cuts", "[oracle]") {
    auto prob = make_small_problem();
    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();

    cptp::sep::SeparationOracle oracle(prob);
    oracle.add_default_separators();

    // All zeros — nothing selected, nothing to violate
    std::vector<double> x_values(m, 0.0);
    std::vector<double> y_values(n, 0.0);

    auto cuts = oracle.separate(x_values, y_values, 0, m, 1e-6);
    REQUIRE(cuts.empty());
}

// =====================================================================
// BoundPropagator: additional tests
// =====================================================================

TEST_CASE("BoundPropagator: has_all_pairs_bounds initially false", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);
    REQUIRE_FALSE(prop.has_all_pairs_bounds());
}

TEST_CASE("BoundPropagator: sweep with loose UB fixes nothing", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    std::vector<double> col_upper = {1.0, 1.0, 1.0};
    // Very loose UB — nothing should be fixed
    double upper_bound = 1e12;

    auto fixings = prop.sweep(upper_bound, col_upper);
    REQUIRE(fixings.empty());
}

TEST_CASE("BoundPropagator: sweep skips already-fixed edges", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 100};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    // Edge 2 is already fixed to 0
    std::vector<double> col_upper = {1.0, 1.0, 0.0};
    double upper_bound = 4.0;

    auto fixings = prop.sweep(upper_bound, col_upper);
    // Edge 2 should NOT appear in fixings because it's already fixed
    for (int32_t e : fixings) {
        REQUIRE(e != 2);
    }
}

TEST_CASE("BoundPropagator: sweep_nodes no fixings when edges are open", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 0, 0};
    std::vector<double> demands = {0, 0, 0};
    prob.build(3, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    int32_t m = prob.num_edges();
    int32_t n = prob.num_nodes();
    // All edges and nodes open
    std::vector<double> col_upper(m + n, 1.0);

    auto node_fixings = prop.sweep_nodes(col_upper, m);
    // All nodes have incident edges open, so no node should be fixable
    REQUIRE(node_fixings.empty());
}

TEST_CASE("BoundPropagator: propagate_fixed_edge does nothing with loose UB", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);

    cptp::preprocess::BoundPropagator prop(prob, fwd, fwd, 0.0);

    int32_t m = prob.num_edges();
    std::vector<double> col_upper(m, 1.0);

    // With a very loose UB, fixing one edge should not cascade
    double upper_bound = 1e12;
    auto fixings = prop.propagate_fixed_edge(0, upper_bound, col_upper);
    REQUIRE(fixings.empty());
}

TEST_CASE("BoundPropagator: accessors return correct values", "[propagator]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {1, 1, 1};
    std::vector<double> profits = {0, 5, 3};
    std::vector<double> demands = {0, 2, 1};
    prob.build(3, edges, costs, profits, demands, 10.0, 0, 0);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);
    auto bwd = cptp::preprocess::backward_labeling(prob, 0);
    double correction = 42.0;

    cptp::preprocess::BoundPropagator prop(prob, fwd, bwd, correction);

    REQUIRE(prop.fwd_bounds().size() == 3);
    REQUIRE(prop.bwd_bounds().size() == 3);
    REQUIRE(prop.correction() == 42.0);
    REQUIRE_FALSE(prop.has_all_pairs_bounds());
}

TEST_CASE("BoundPropagator: sweep on s-t path problem", "[propagator][path]") {
    // 4-node path, source=0, target=3. Edge {1,2} cost=100.
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {1, 1, 1, 100, 1, 1};
    std::vector<double> profits = {0, 0, 0, 0};
    std::vector<double> demands = {0, 0, 0, 0};
    prob.build(4, edges, costs, profits, demands, 1e18, 0, 3);

    auto fwd = cptp::preprocess::forward_labeling(prob, 0);
    auto bwd = cptp::preprocess::backward_labeling(prob, 3);

    cptp::preprocess::BoundPropagator prop(prob, fwd, bwd, 0.0);

    int32_t m = prob.num_edges();
    std::vector<double> col_upper(m, 1.0);
    double upper_bound = 5.0;

    auto fixings = prop.sweep(upper_bound, col_upper);

    // Edge {1,2} (index 3) with cost=100 should be fixed
    bool fixes_expensive = false;
    for (int32_t e : fixings) {
        if (e == 3) fixes_expensive = true;
    }
    REQUIRE(fixes_expensive);
}
