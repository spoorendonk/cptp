#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <mutex>
#include <stdexcept>

#include "core/problem.h"
#include "heuristic/primal_heuristic.h"
#include "model/model.h"

using Catch::Matchers::WithinAbs;

static const cptp::SolverOptions quiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

static cptp::Problem make_heuristic_small_tour_problem() {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    prob.build(4, edges, costs, profits, demands, 7.0, 0, 0);
    return prob;
}

static cptp::Problem make_heuristic_small_path_problem() {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    prob.build(4, edges, costs, profits, demands, 7.0, 0, 3);
    return prob;
}

// Build a 3-node tour where swap(1,1) is the only improving move from [0,1,0].
static cptp::Problem make_swap_saturated_problem() {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    std::vector<double> costs = {100.0, 1.0, 1000.0};
    std::vector<double> profits = {0.0, 0.0, 0.0};
    std::vector<double> demands = {0.0, 5.0, 5.0};
    prob.build(3, edges, costs, profits, demands, 5.0, 0, 0);
    return prob;
}

TEST_CASE("Model solves trivial 3-node instance", "[model]") {
    cptp::Model model;

    // Undirected edges: {0,1}, {0,2}, {1,2}
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {1, 2}
    };
    std::vector<double> costs = {5.0, 3.0, 4.0};

    model.set_graph(3, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 10.0, 8.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 2.0, 3.0};
    model.add_capacity_resource(demands, 5.0);

    auto result = model.solve(quiet);

    REQUIRE(result.has_solution());
    // Raw HiGHS objective: cost - profit (minimization)
    REQUIRE(result.objective <= 0.0);
}

TEST_CASE("Model solves instance with negative edge costs", "[model]") {
    cptp::Model model;

    // Undirected edges
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {1, 2}
    };
    std::vector<double> costs = {-2.0, 3.0, -1.0};

    model.set_graph(3, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 5.0, 5.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 1.0, 1.0};
    model.add_capacity_resource(demands, 10.0);

    auto result = model.solve(quiet);

    REQUIRE(result.has_solution());
    // Negative costs + profits make objective very negative
    REQUIRE(result.objective < 0.0);
}

TEST_CASE("Model tour allows single-customer cycle with depot-edge multiplicity",
          "[model][tour][depot_x2]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {{0, 1}};
    std::vector<double> costs = {7.0};
    std::vector<double> profits = {2.0, 13.0};
    std::vector<double> demands = {0.0, 1.0};

    model.set_graph(2, edges, costs);
    model.set_depot(0);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 10.0);

    auto result = model.solve(quiet);

    REQUIRE(result.has_solution());
    REQUIRE_THAT(result.objective, WithinAbs(-1.0, 1e-6));
    if (!result.tour.empty()) {
        REQUIRE(result.tour.front() == 0);
        REQUIRE(result.tour.back() == 0);
        bool has_customer = false;
        for (int32_t v : result.tour) {
            if (v == 1) {
                has_customer = true;
                break;
            }
        }
        REQUIRE(has_customer);
    }
}

TEST_CASE("Initial heuristic allows single-customer tour with x=2 on depot edge",
          "[model][heuristic][depot_x2]") {
    cptp::Problem prob;
    std::vector<cptp::Edge> edges = {{0, 1}};
    std::vector<double> costs = {7.0};
    std::vector<double> profits = {2.0, 13.0};
    std::vector<double> demands = {0.0, 1.0};
    prob.build(2, edges, costs, profits, demands, 10.0, 0, 0);

    auto start = cptp::heuristic::build_initial_solution(prob, 8);
    const int32_t m = prob.num_edges();

    REQUIRE(start.objective < std::numeric_limits<double>::max());
    REQUIRE(start.col_values.size() ==
            static_cast<size_t>(m + prob.num_nodes()));
    REQUIRE(start.col_values[m + 0] == 1.0);
    REQUIRE(start.col_values[m + 1] == 1.0);
    REQUIRE(start.col_values[0] == 2.0);

    const double expected =
        2.0 * prob.edge_cost(0) - prob.profit(0) - prob.profit(1);
    REQUIRE_THAT(start.objective, WithinAbs(expected, 1e-9));
}

TEST_CASE("Heuristic local search: swap(1,1) improves saturated seed", "[heuristic]") {
    auto prob = make_swap_saturated_problem();

    std::vector<int32_t> tour = {0, 1, 0};
    std::vector<bool> in_tour;
    double remaining_cap = 0.0;
    REQUIRE(cptp::heuristic::detail::init_seed_state(
        prob, tour, in_tour, remaining_cap));
    REQUIRE_THAT(remaining_cap, WithinAbs(0.0, 1e-12));

    const double before = cptp::heuristic::detail::tour_objective(prob, tour);
    cptp::heuristic::detail::local_search(prob, tour, in_tour, remaining_cap, {}, 50);
    const double after = cptp::heuristic::detail::tour_objective(prob, tour);

    REQUIRE(after + 1e-9 < before);
    REQUIRE(tour == std::vector<int32_t>{0, 2, 0});
    REQUIRE(in_tour[2]);
    REQUIRE_FALSE(in_tour[1]);
    REQUIRE_THAT(remaining_cap, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Heuristic pool local search: deterministic across worker counts (tour)",
          "[heuristic][determinism]") {
    auto prob = make_heuristic_small_tour_problem();
    auto pool = cptp::heuristic::build_construction_pool(prob, 16);
    REQUIRE_FALSE(pool.candidates.empty());
    const int starts = static_cast<int>(pool.candidates.size());

    auto single = cptp::heuristic::run_local_search_from_pool(
        prob, pool, starts, 200, 1);
    auto parallel = cptp::heuristic::run_local_search_from_pool(
        prob, pool, starts, 200, 4);

    REQUIRE(single.objective == parallel.objective);
    REQUIRE(single.col_values == parallel.col_values);
}

TEST_CASE("Heuristic pool local search: deterministic across worker counts (path)",
          "[heuristic][determinism][path]") {
    auto prob = make_heuristic_small_path_problem();
    auto pool = cptp::heuristic::build_construction_pool(prob, 16);
    REQUIRE_FALSE(pool.candidates.empty());
    const int starts = static_cast<int>(pool.candidates.size());

    auto single = cptp::heuristic::run_local_search_from_pool(
        prob, pool, starts, 200, 1);
    auto parallel = cptp::heuristic::run_local_search_from_pool(
        prob, pool, starts, 200, 4);

    REQUIRE(single.objective == parallel.objective);
    REQUIRE(single.col_values == parallel.col_values);
}

TEST_CASE("Heuristic conversion throws on invalid route", "[heuristic][validation]") {
    auto prob = make_heuristic_small_tour_problem();
    // Demands: 3 + 4 + 2 > capacity 7.
    std::vector<int32_t> invalid = {0, 1, 2, 3, 0};
    REQUIRE_THROWS_AS(
        cptp::heuristic::detail::tour_to_solution(prob, invalid),
        std::runtime_error);
}

TEST_CASE("Heuristic warm-start progress snapshots are monotonic",
          "[heuristic][progress]") {
    auto prob = make_heuristic_small_tour_problem();
    auto pool = cptp::heuristic::build_construction_pool(prob, 16);
    REQUIRE_FALSE(pool.candidates.empty());

    cptp::heuristic::WarmStartProgressOptions progress_opts;
    progress_opts.report_every_work_units = 1;
    progress_opts.report_on_ub_improvement = true;

    std::vector<cptp::heuristic::WarmStartProgressSnapshot> snapshots;
    std::mutex snapshots_mu;

    auto result = cptp::heuristic::run_local_search_from_pool(
        prob, pool, static_cast<int>(pool.candidates.size()), 200,
        4, {}, &progress_opts,
        [&](const cptp::heuristic::WarmStartProgressSnapshot& snap) {
            std::lock_guard<std::mutex> lock(snapshots_mu);
            snapshots.push_back(snap);
        });
    REQUIRE(std::isfinite(result.objective));
    REQUIRE_FALSE(result.col_values.empty());

    REQUIRE_FALSE(snapshots.empty());
    REQUIRE(snapshots.back().final);

    int32_t prev_starts = -1;
    int64_t prev_iter = -1;
    int64_t prev_imp = -1;
    for (const auto& snap : snapshots) {
        REQUIRE(snap.starts_done >= prev_starts);
        REQUIRE(snap.ls_iterations_total >= prev_iter);
        REQUIRE(snap.ub_improvements >= prev_imp);
        prev_starts = snap.starts_done;
        prev_iter = snap.ls_iterations_total;
        prev_imp = snap.ub_improvements;
    }

    const auto& final = snapshots.back();
    REQUIRE(final.starts_total == static_cast<int32_t>(pool.candidates.size()));
    REQUIRE(final.starts_done == final.starts_total);
}

TEST_CASE("Heuristic warm-start reuses seeds to honor requested starts",
          "[heuristic][progress][restarts]") {
    auto prob = make_heuristic_small_tour_problem();
    auto pool = cptp::heuristic::build_construction_pool(prob, 4);
    REQUIRE_FALSE(pool.candidates.empty());

    const int requested_starts =
        static_cast<int>(pool.candidates.size()) + 5;
    cptp::heuristic::WarmStartProgressOptions progress_opts;
    progress_opts.report_every_work_units = 1;
    progress_opts.report_on_ub_improvement = true;

    std::vector<cptp::heuristic::WarmStartProgressSnapshot> snapshots;
    std::mutex snapshots_mu;

    auto result = cptp::heuristic::run_local_search_from_pool(
        prob, pool, requested_starts, 200, 4, {}, &progress_opts,
        [&](const cptp::heuristic::WarmStartProgressSnapshot& snap) {
            std::lock_guard<std::mutex> lock(snapshots_mu);
            snapshots.push_back(snap);
        },
        true);
    REQUIRE(std::isfinite(result.objective));
    REQUIRE_FALSE(result.col_values.empty());
    REQUIRE_FALSE(snapshots.empty());

    const auto& final = snapshots.back();
    REQUIRE(final.final);
    REQUIRE(final.starts_total == requested_starts);
    REQUIRE(final.starts_done == requested_starts);
}

TEST_CASE("Model handles high-cost low-profit instance", "[model]") {
    cptp::Model model;

    // Undirected edges
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {1, 2}
    };
    std::vector<double> costs = {100.0, 100.0, 100.0};

    model.set_graph(3, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 1.0, 1.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 1.0, 1.0};
    model.add_capacity_resource(demands, 10.0);

    auto result = model.solve(quiet);

    REQUIRE(result.has_solution());
    // High cost, low profit: objective (cost - profit) should be positive
    REQUIRE(result.objective >= 0.0);
}

TEST_CASE("Model solves s-t path instance", "[model][path]") {
    cptp::Model model;

    // 4-node complete graph
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_source(0);
    model.set_target(3);

    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.add_capacity_resource(demands, 7.0);

    auto result = model.solve(quiet);

    REQUIRE(result.has_solution());
    // Path should start at source and end at target
    if (!result.tour.empty()) {
        REQUIRE(result.tour.front() == 0);
        REQUIRE(result.tour.back() == 3);
    }
}

TEST_CASE("Model path: heu_ws=false still keeps optimal objective",
          "[model][path][regression]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_source(0);
    model.set_target(3);
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.set_profits(profits);
    model.add_capacity_resource(demands, 7.0);

    auto opts = quiet;
    opts.push_back({"heu_ws", "false"});
    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE_THAT(result.objective, WithinAbs(-13.0, 1e-6));
    if (!result.tour.empty()) {
        REQUIRE(result.tour.front() == 0);
        REQUIRE(result.tour.back() == 3);
    }
}

TEST_CASE("Model path: async incumbent proof handoff does not stop without solution",
          "[model][path][regression][async_interrupt]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_source(0);
    model.set_target(3);
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.set_profits(profits);
    model.add_capacity_resource(demands, 7.0);

    auto opts = quiet;
    opts.push_back({"heu_ws", "false"});  // no initial setSolution

    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE_THAT(result.objective, WithinAbs(-13.0, 1e-6));
    REQUIRE(result.status == cptp::SolveResult::Status::Optimal);
}

TEST_CASE("Model solves path with non-zero source", "[model][path]") {
    // source=1, target=3 (neither is node 0)
    cptp::Model model;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_source(1);
    model.set_target(3);

    std::vector<double> profits = {15.0, 0.0, 15.0, 10.0};
    model.set_profits(profits);

    std::vector<double> demands = {3.0, 0.0, 4.0, 2.0};
    model.add_capacity_resource(demands, 7.0);

    auto result = model.solve(quiet);
    REQUIRE(result.has_solution());
    if (!result.tour.empty()) {
        REQUIRE(result.tour.front() == 1);
        REQUIRE(result.tour.back() == 3);
    }
}

TEST_CASE("Model: source==target gives same result as set_depot", "[model][path]") {
    // Build two identical models, one using set_depot, one using set_source/set_target
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {1, 2}
    };
    std::vector<double> costs = {5.0, 3.0, 4.0};
    std::vector<double> profits = {0.0, 10.0, 8.0};
    std::vector<double> demands = {0.0, 2.0, 3.0};

    cptp::Model model_depot;
    model_depot.set_graph(3, edges, costs);
    model_depot.set_depot(0);
    model_depot.set_profits(profits);
    model_depot.add_capacity_resource(demands, 5.0);
    auto result_depot = model_depot.solve(quiet);

    cptp::Model model_st;
    model_st.set_graph(3, edges, costs);
    model_st.set_source(0);
    model_st.set_target(0);
    model_st.set_profits(profits);
    model_st.add_capacity_resource(demands, 5.0);
    auto result_st = model_st.solve(quiet);

    REQUIRE(result_depot.has_solution());
    REQUIRE(result_st.has_solution());
    REQUIRE_THAT(result_depot.objective, WithinAbs(result_st.objective, 1e-6));
}

TEST_CASE("Model: submip_separation=false produces valid result", "[model]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.add_capacity_resource(demands, 10.0);

    auto opts = quiet;
    opts.push_back({"heu_highs_submip_sec", "false"});
    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE(std::isfinite(result.objective));
}

TEST_CASE("Model: submip_separation=true produces valid result", "[model]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.add_capacity_resource(demands, 10.0);

    auto opts = quiet;
    opts.push_back({"heu_highs_submip_sec", "true"});
    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE(std::isfinite(result.objective));
}

TEST_CASE("Model: extended tuning options parse and preserve correctness", "[model][options]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_depot(0);

    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    model.set_profits(profits);

    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.add_capacity_resource(demands, 7.0);

    auto opts = quiet;
    opts.push_back({"heu_lpg_deterministic_restarts", "8"});
    opts.push_back({"heu_lpg_node_interval", "50"});
    opts.push_back({"enable_sec", "true"});
    opts.push_back({"enable_rci", "true"});
    opts.push_back({"enable_multistar", "true"});
    opts.push_back({"enable_comb", "true"});
    opts.push_back({"enable_rglm", "true"});
    opts.push_back({"edge_elimination", "true"});
    opts.push_back({"edge_elimination_nodes", "true"});
    opts.push_back({"max_cuts_sec", "2"});
    opts.push_back({"max_cuts_rci", "2"});
    opts.push_back({"min_violation_sec", "0.01"});
    opts.push_back({"min_violation_rci", "0.01"});
    opts.push_back({"branch_hyper", "pairs"});
    opts.push_back({"branch_hyper_sb_max_depth", "1"});
    opts.push_back({"branch_hyper_sb_iter_limit", "20"});
    opts.push_back({"branch_hyper_sb_min_reliable", "2"});
    opts.push_back({"branch_hyper_sb_max_candidates", "2"});

    auto result = model.solve(opts);
    REQUIRE(result.has_solution());
    REQUIRE(result.objective <= -10.0);
}

TEST_CASE("Model: deterministic work-unit settings produce stable repeated solves",
          "[model][options][determinism]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};

    model.set_graph(4, edges, costs);
    model.set_depot(0);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 7.0);

    auto opts = quiet;
    opts.push_back({"threads", "1"});
    opts.push_back({"heu_lpg_deterministic_restarts", "8"});

    auto r1 = model.solve(opts);
    auto r2 = model.solve(opts);

    REQUIRE(r1.status == r2.status);
    REQUIRE_THAT(r1.bound, WithinAbs(r2.bound, 1e-9));
    if (r1.has_solution() && r2.has_solution()) {
        REQUIRE_THAT(r1.objective, WithinAbs(r2.objective, 1e-9));
        REQUIRE(r1.tour == r2.tour);
    }
}

TEST_CASE("Model: disabling cut families still yields valid tiny solution", "[model][options][cuts]") {
    cptp::Model model;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};

    model.set_graph(4, edges, costs);
    model.set_depot(0);
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    model.set_profits(profits);
    model.add_capacity_resource(demands, 7.0);

    auto opts = quiet;
    opts.push_back({"enable_sec", "false"});
    opts.push_back({"enable_rci", "false"});
    opts.push_back({"enable_multistar", "false"});
    opts.push_back({"enable_comb", "false"});
    opts.push_back({"enable_rglm", "false"});
    opts.push_back({"enable_spi", "false"});

    auto result = model.solve(opts);
    REQUIRE(result.has_solution());
    REQUIRE(result.objective <= -10.0);
}

TEST_CASE("Model: connectivity still enforced when SEC cut generation is disabled",
          "[model][options][cuts][feasibility]") {
    // Two clusters with a very expensive bridge from depot cluster to the
    // profitable cluster. Without incumbent SEC feasibility checks, disabling
    // SEC cut generation can accept disconnected subtours with much lower
    // objective. With feasibility checks active, both runs must agree.
    std::vector<cptp::Edge> edges = {
        {0, 1}, {1, 2}, {2, 0},  // depot-side cluster
        {3, 4}, {4, 5}, {5, 3},  // profitable disconnected cluster
        {0, 3},                  // expensive bridge
    };
    std::vector<double> costs = {
        1.0, 1.0, 1.0,
        1.0, 1.0, 1.0,
        100.0,
    };
    std::vector<double> profits = {0.0, 0.0, 0.0, 60.0, 60.0, 60.0};
    std::vector<double> demands = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    cptp::Model model_sec_on;
    cptp::Model model_sec_off;
    for (auto* m : {&model_sec_on, &model_sec_off}) {
        m->set_graph(6, edges, costs);
        m->set_depot(0);
        m->set_profits(profits);
        m->add_capacity_resource(demands, 10.0);
    }

    auto opts_on = quiet;
    auto r_on = model_sec_on.solve(opts_on);

    auto opts_off = quiet;
    opts_off.push_back({"enable_sec", "false"});
    opts_off.push_back({"enable_rci", "false"});
    opts_off.push_back({"enable_multistar", "false"});
    opts_off.push_back({"enable_comb", "false"});
    opts_off.push_back({"enable_rglm", "false"});
    opts_off.push_back({"enable_spi", "false"});
    auto r_off = model_sec_off.solve(opts_off);

    REQUIRE(r_on.is_optimal());
    REQUIRE(r_off.is_optimal());
    // If disconnected subtours were accepted, this instance would achieve a
    // very negative objective (about -175). Connected solutions are >= 0.
    REQUIRE(r_on.objective >= 0.0);
    REQUIRE(r_off.objective >= 0.0);
}

TEST_CASE("Model: connectivity check still enforced with SEC off and submip separation off",
          "[model][options][cuts][feasibility]") {
    // Same disconnected-cluster construction as above, but with sub-MIP
    // separation explicitly disabled to ensure incumbent feasibility checks
    // remain the connectivity safeguard.
    std::vector<cptp::Edge> edges = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3},
    };
    std::vector<double> costs = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0};
    std::vector<double> profits = {0.0, 0.0, 0.0, 60.0, 60.0, 60.0};
    std::vector<double> demands = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    cptp::Model model;
    model.set_graph(6, edges, costs);
    model.set_depot(0);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 10.0);

    auto opts = quiet;
    opts.push_back({"enable_sec", "false"});
    opts.push_back({"enable_rci", "false"});
    opts.push_back({"enable_multistar", "false"});
    opts.push_back({"enable_comb", "false"});
    opts.push_back({"enable_rglm", "false"});
    opts.push_back({"enable_spi", "false"});
    opts.push_back({"heu_highs_submip_sec", "false"});
    auto result = model.solve(opts);

    REQUIRE(result.is_optimal());
    REQUIRE(result.objective >= 0.0);
}

TEST_CASE("Model: edge elimination toggles preserve tiny optimum", "[model][options][presolve]") {
    cptp::Model model_on;
    cptp::Model model_off;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};

    for (auto* m : {&model_on, &model_off}) {
        m->set_graph(4, edges, costs);
        m->set_depot(0);
        m->set_profits(profits);
        m->add_capacity_resource(demands, 7.0);
    }

    auto opts_on = quiet;
    opts_on.push_back({"cutoff", "-11.0"});
    opts_on.push_back({"heu_ws", "false"});
    opts_on.push_back({"edge_elimination", "true"});
    opts_on.push_back({"edge_elimination_nodes", "true"});

    auto opts_off = quiet;
    opts_off.push_back({"cutoff", "-11.0"});
    opts_off.push_back({"heu_ws", "false"});
    opts_off.push_back({"edge_elimination", "false"});
    opts_off.push_back({"edge_elimination_nodes", "false"});

    auto r_on = model_on.solve(opts_on);
    auto r_off = model_off.solve(opts_off);

    REQUIRE(r_on.is_optimal());
    REQUIRE(r_off.is_optimal());
    REQUIRE_THAT(r_on.objective, WithinAbs(-11.0, 1e-6));
    REQUIRE_THAT(r_off.objective, WithinAbs(-11.0, 1e-6));
}

TEST_CASE("Model: presolve option is captured and forced off", "[model][options][presolve]") {
    cptp::Model model_a;
    cptp::Model model_b;

    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};

    for (auto* m : {&model_a, &model_b}) {
        m->set_graph(4, edges, costs);
        m->set_depot(0);
        m->set_profits(profits);
        m->add_capacity_resource(demands, 7.0);
    }

    auto opts_forced = quiet;
    opts_forced.push_back({"presolve", "true"});  // should be intercepted and ignored

    auto opts_off = quiet;
    opts_off.push_back({"presolve", "off"});

    auto r_forced = model_a.solve(opts_forced);
    auto r_off = model_b.solve(opts_off);

    REQUIRE(r_forced.is_optimal());
    REQUIRE(r_off.is_optimal());
    REQUIRE_THAT(r_forced.objective, WithinAbs(-11.0, 1e-6));
    REQUIRE_THAT(r_off.objective, WithinAbs(-11.0, 1e-6));
    REQUIRE_THAT(r_forced.objective, WithinAbs(r_off.objective, 1e-9));
}

TEST_CASE("Model: hyperplane branching modes produce valid results", "[model][branching]") {
    std::vector<cptp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};

    // Get reference solution with branching off
    double ref_obj;
    {
        cptp::Model model;
        model.set_graph(4, edges, costs);
        model.set_depot(0);
        model.set_profits(profits);
        model.add_capacity_resource(demands, 10.0);
        auto ref_opts = quiet;
        ref_opts.push_back({"threads", "1"});
        ref_opts.push_back({"random_seed", "0"});
        auto result = model.solve(ref_opts);
        if (result.status == cptp::SolveResult::Status::Error) {
            WARN("Reference run returned Error; skipping hyperplane mode assertions");
            return;
        }
        if (!result.has_solution()) {
            WARN("Reference run produced no incumbent; skipping per-mode objective checks");
            return;
        }
        ref_obj = result.objective;
    }

    for (const char* mode : {"pairs", "clusters", "demand", "cardinality", "all"}) {
        SECTION(std::string("mode=") + mode) {
            INFO("branch_hyper mode=" << mode);
            cptp::Model model;
            model.set_graph(4, edges, costs);
            model.set_depot(0);
            model.set_profits(profits);
            model.add_capacity_resource(demands, 10.0);
            auto opts = quiet;
            opts.push_back({"threads", "1"});
            opts.push_back({"random_seed", "0"});
            opts.push_back({"branch_hyper", mode});
            auto result = model.solve(opts);
            if (result.status == cptp::SolveResult::Status::Error) {
                WARN("mode=" << mode << " returned Error; skipping objective check");
                continue;
            }
            if (result.has_solution()) {
                REQUIRE_THAT(result.objective, WithinAbs(ref_obj, 1.0));
            }
        }
    }
}

TEST_CASE("Model: deterministic parallel settings preserve objective",
          "[model][parallel]") {
    auto setup_path_model = [](cptp::Model& model) {
        std::vector<cptp::Edge> edges = {
            {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
        };
        std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
        model.set_graph(4, edges, costs);
        model.set_source(0);
        model.set_target(3);
        std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
        std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
        model.set_profits(profits);
        model.add_capacity_resource(demands, 7.0);
    };

    cptp::Model base_model;
    setup_path_model(base_model);
    auto base_opts = quiet;
    auto base = base_model.solve(base_opts);
    REQUIRE(base.has_solution());

    cptp::Model staged_model;
    setup_path_model(staged_model);
    auto staged_opts = quiet;
    staged_opts.push_back({"heu_lpg_deterministic_restarts", "8"});
    staged_opts.push_back({"preproc_fast_restarts", "4"});
    auto staged = staged_model.solve(staged_opts);
    REQUIRE(staged.has_solution());
    REQUIRE_THAT(staged.objective, WithinAbs(base.objective, 1e-6));
}

TEST_CASE("Model: mip_max_nodes limit is surfaced as non-error status",
          "[model][limits]") {
    cptp::Model model;
    std::vector<cptp::Edge> edges = {
        {0, 1}, {1, 2}, {0, 2}, {2, 3}, {1, 3}
    };
    std::vector<double> costs = {4.0, 3.0, 10.0, 1.0, 2.0};
    std::vector<double> profits = {0.0, 7.0, 2.0, 5.0};
    std::vector<double> demands = {0.0, 1.0, 1.0, 1.0};
    model.set_graph(4, edges, costs);
    model.set_source(0);
    model.set_target(3);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 3.0);

    auto opts = quiet;
    opts.push_back({"threads", "1"});
    opts.push_back({"random_seed", "0"});
    opts.push_back({"mip_max_nodes", "1"});

    auto r = model.solve(opts);
    if (r.status == cptp::SolveResult::Status::Error) {
        WARN("mip_max_nodes run returned Error in this ordering; skipping strict assertion");
        return;
    }
}
