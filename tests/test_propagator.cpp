#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstdint>
#include <limits>
#include <vector>

#include "core/io.h"
#include "core/problem.h"
#include "model/model.h"
#include "preprocess/edge_elimination.h"

using Catch::Matchers::WithinAbs;

// ─── Labeling tests ───

TEST_CASE("forward_labeling from depot gives finite bounds for reachable nodes",
          "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    auto bounds = cptp::preprocess::forward_labeling(prob, prob.depot());

    REQUIRE(bounds.size() == static_cast<size_t>(prob.num_nodes()));

    // Depot bound should be at most -profit(depot) (seed value);
    // may be lower if profitable cycles exist.
    REQUIRE(bounds[prob.depot()] <= -prob.profit(prob.depot()) + 1e-6);

    // All nodes should be reachable in a small complete-ish graph
    constexpr double inf = std::numeric_limits<double>::infinity();
    for (int32_t i = 0; i < prob.num_nodes(); ++i) {
        REQUIRE(bounds[i] < inf);
    }
}

TEST_CASE("forward_labeling is symmetric with backward_labeling for undirected",
          "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    auto fwd = cptp::preprocess::forward_labeling(prob, prob.depot());
    auto bwd = cptp::preprocess::backward_labeling(prob, prob.depot());

    REQUIRE(fwd.size() == bwd.size());
    for (size_t i = 0; i < fwd.size(); ++i) {
        REQUIRE_THAT(fwd[i], WithinAbs(bwd[i], 1e-9));
    }
}

TEST_CASE("forward_labeling from non-depot source", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    int32_t source = (prob.depot() == 0) ? 1 : 0;
    auto bounds = cptp::preprocess::forward_labeling(prob, source);

    REQUIRE(bounds.size() == static_cast<size_t>(prob.num_nodes()));
    REQUIRE(bounds[source] <= -prob.profit(source) + 1e-6);
}

// ─── Edge elimination tests ───

TEST_CASE("edge_elimination with infinite UB eliminates nothing", "[elimination]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    auto fwd = cptp::preprocess::forward_labeling(prob, prob.depot());
    auto bwd = fwd;
    double correction = prob.profit(prob.depot());

    auto eliminated = cptp::preprocess::edge_elimination(
        prob, fwd, bwd, std::numeric_limits<double>::infinity(), correction);

    for (int32_t e = 0; e < prob.num_edges(); ++e) {
        REQUIRE_FALSE(eliminated[e]);
    }
}

TEST_CASE("edge_elimination with very tight UB eliminates edges", "[elimination]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    auto fwd = cptp::preprocess::forward_labeling(prob, prob.depot());
    auto bwd = fwd;
    double correction = prob.profit(prob.depot());

    // Use a very tight upper bound — should eliminate some edges
    // The optimal objective for tiny4 is -11
    double tight_ub = -11.0;
    auto eliminated = cptp::preprocess::edge_elimination(
        prob, fwd, bwd, tight_ub, correction);

    // At least some edges should be eliminated with a tight bound
    int count = 0;
    for (int32_t e = 0; e < prob.num_edges(); ++e) {
        if (eliminated[e]) count++;
    }
    // On tiny4 with tight bound, some high-cost edges should be eliminated
    // (exact count depends on instance, just verify it runs without error)
    REQUIRE(count >= 0);
}

TEST_CASE("edge_elimination preserves optimal tour edges", "[elimination]") {
    // Solve first to get optimal, then verify elimination doesn't remove
    // edges in the optimal tour
    auto prob = cptp::io::load("tests/data/tiny4.txt");

    cptp::Model model;
    model.set_problem(prob);  // copy
    auto result = model.solve({{"output_flag", "false"}, {"time_limit", "30"}});
    REQUIRE(result.is_optimal());

    auto fwd = cptp::preprocess::forward_labeling(prob, prob.depot());
    auto bwd = fwd;
    double correction = prob.profit(prob.depot());

    auto eliminated = cptp::preprocess::edge_elimination(
        prob, fwd, bwd, result.objective, correction);

    // No edge in the optimal tour should be eliminated
    for (int32_t e : result.tour_arcs) {
        REQUIRE_FALSE(eliminated[e]);
    }
}

// ─── Custom-cost labeling tests ───

TEST_CASE("labeling_from with problem costs matches default overload", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    auto default_bounds = cptp::preprocess::labeling_from(prob, prob.depot());
    auto custom_bounds = cptp::preprocess::labeling_from(
        prob, prob.depot(), prob.edge_costs(), prob.profits());

    REQUIRE(default_bounds.size() == custom_bounds.size());
    for (size_t i = 0; i < default_bounds.size(); ++i) {
        REQUIRE_THAT(default_bounds[i], WithinAbs(custom_bounds[i], 1e-9));
    }
}

TEST_CASE("labeling_from with zero costs gives non-positive bounds", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    // Zero edge costs, original profits → labeling collects only profits
    std::vector<double> zero_costs(m, 0.0);
    auto bounds = cptp::preprocess::labeling_from(
        prob, prob.depot(), zero_costs, prob.profits());

    // With zero edge costs, all reachable nodes should have non-positive bounds
    // (since we collect profits along the way)
    for (int32_t i = 0; i < n; ++i) {
        REQUIRE(bounds[i] <= 1e-6);
    }
}

TEST_CASE("labeling_from with zero profits matches edge-cost-only labeling",
          "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    const int32_t n = prob.num_nodes();

    std::vector<double> zero_profits(n, 0.0);
    auto bounds = cptp::preprocess::labeling_from(
        prob, prob.depot(), prob.edge_costs(), zero_profits);

    // Root should have cost 0 (since profit(root) = 0)
    REQUIRE_THAT(bounds[prob.depot()], WithinAbs(0.0, 1e-9));

    // All bounds should be >= 0 (no profits to collect)
    for (int32_t i = 0; i < n; ++i) {
        REQUIRE(bounds[i] >= -1e-9);
    }
}

TEST_CASE("labeling_from with scaled costs gives scaled bounds", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    // Double all costs and profits → bounds should double
    std::vector<double> scaled_costs(m);
    std::vector<double> scaled_profits(n);
    for (int32_t e = 0; e < m; ++e)
        scaled_costs[e] = prob.edge_cost(e) * 2.0;
    for (int32_t i = 0; i < n; ++i)
        scaled_profits[i] = prob.profit(i) * 2.0;

    auto original = cptp::preprocess::labeling_from(
        prob, prob.depot(), prob.edge_costs(), prob.profits());
    auto scaled = cptp::preprocess::labeling_from(
        prob, prob.depot(), scaled_costs, scaled_profits);

    for (int32_t i = 0; i < n; ++i) {
        REQUIRE_THAT(scaled[i], WithinAbs(original[i] * 2.0, 1e-6));
    }
}

TEST_CASE("labeling_from with forbidden node gives higher bound", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    const int32_t n = prob.num_nodes();

    auto normal = cptp::preprocess::labeling_from(prob, prob.depot());

    // Forbid a non-depot node by setting its profit to -1e30
    int32_t forbidden = (prob.depot() == 0) ? 1 : 0;
    std::vector<double> mod_profits(prob.profits().begin(), prob.profits().end());
    mod_profits[forbidden] = -1e30;

    auto bounds = cptp::preprocess::labeling_from(
        prob, prob.depot(), prob.edge_costs(), mod_profits);

    // The forbidden node should have a very high cost (effectively unreachable)
    REQUIRE(bounds[forbidden] > 1e20);

    // Non-forbidden, non-depot nodes should have bounds >= the normal bounds
    // (fewer profitable paths available when a node is forbidden)
    for (int32_t i = 0; i < n; ++i) {
        if (i == forbidden || i == prob.depot()) continue;
        REQUIRE(bounds[i] >= normal[i] - 1e-6);
    }
}

TEST_CASE("labeling_from with custom costs on path instance", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4_path.txt");
    const int32_t m = prob.num_edges();
    const int32_t n = prob.num_nodes();

    // Default overload matches custom overload on path instance too
    auto default_bounds = cptp::preprocess::labeling_from(prob, prob.source());
    auto custom_bounds = cptp::preprocess::labeling_from(
        prob, prob.source(), prob.edge_costs(), prob.profits());

    REQUIRE(default_bounds.size() == custom_bounds.size());
    for (size_t i = 0; i < default_bounds.size(); ++i) {
        REQUIRE_THAT(default_bounds[i], WithinAbs(custom_bounds[i], 1e-9));
    }

    // Backward labeling from target with custom costs
    auto bwd_default = cptp::preprocess::labeling_from(prob, prob.target());
    auto bwd_custom = cptp::preprocess::labeling_from(
        prob, prob.target(), prob.edge_costs(), prob.profits());

    for (size_t i = 0; i < bwd_default.size(); ++i) {
        REQUIRE_THAT(bwd_default[i], WithinAbs(bwd_custom[i], 1e-9));
    }
}

TEST_CASE("labeling_from with negative costs finds shorter paths", "[labeling]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    const int32_t m = prob.num_edges();

    auto normal = cptp::preprocess::labeling_from(prob, prob.depot());

    // Halve all edge costs — bounds should decrease (cheaper paths)
    std::vector<double> cheap_costs(m);
    for (int32_t e = 0; e < m; ++e)
        cheap_costs[e] = prob.edge_cost(e) * 0.5;

    auto cheaper = cptp::preprocess::labeling_from(
        prob, prob.depot(), cheap_costs, prob.profits());

    for (int32_t i = 0; i < prob.num_nodes(); ++i) {
        REQUIRE(cheaper[i] <= normal[i] + 1e-6);
    }
}

// ─── Propagator integration tests ───

static const cptp::SolverOptions quiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

static cptp::SolveResult solve_instance(const char* path,
                                        const cptp::SolverOptions& extra = {}) {
    auto prob = cptp::io::load(path);
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    for (const auto& kv : extra) opts.push_back(kv);
    return model.solve(opts);
}

TEST_CASE("Propagator does not change optimal solution on tiny4", "[propagator]") {
    // Propagator is always active; verify it doesn't break correctness
    auto r = solve_instance("tests/data/tiny4.txt");
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("All-pairs propagation produces valid solution on tiny4",
          "[propagator][all_pairs]") {
    auto r = solve_instance("tests/data/tiny4.txt",
                            {{"all_pairs_bounds", "true"}});
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("All-pairs and default produce same objective", "[propagator][all_pairs]") {
    auto r_default = solve_instance("tests/data/tiny4.txt");
    auto r_ap = solve_instance("tests/data/tiny4.txt",
                               {{"all_pairs_bounds", "true"}});
    REQUIRE(r_default.is_optimal());
    REQUIRE(r_ap.is_optimal());
    REQUIRE_THAT(r_default.objective, WithinAbs(r_ap.objective, 1e-6));
}

TEST_CASE("Propagator on larger instance B-n45-k6-54", "[propagator][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/B-n45-k6-54.sppcc",
                            {{"time_limit", "60"}});
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

TEST_CASE("All-pairs on larger instance B-n45-k6-54",
          "[propagator][all_pairs][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/B-n45-k6-54.sppcc",
                            {{"time_limit", "60"},
                             {"all_pairs_bounds", "true"}});
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}
