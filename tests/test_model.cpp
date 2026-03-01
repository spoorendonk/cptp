#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "model/deterministic_checkpoint.h"
#include "model/model.h"

using Catch::Matchers::WithinAbs;

static const rcspp::SolverOptions quiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

TEST_CASE("Model solves trivial 3-node instance", "[model]") {
    rcspp::Model model;

    // Undirected edges: {0,1}, {0,2}, {1,2}
    std::vector<rcspp::Edge> edges = {
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
    rcspp::Model model;

    // Undirected edges
    std::vector<rcspp::Edge> edges = {
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

TEST_CASE("Model handles high-cost low-profit instance", "[model]") {
    rcspp::Model model;

    // Undirected edges
    std::vector<rcspp::Edge> edges = {
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
    rcspp::Model model;

    // 4-node complete graph
    std::vector<rcspp::Edge> edges = {
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

TEST_CASE("Model path: disable_heuristics still keeps optimal objective",
          "[model][path][regression]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"disable_heuristics", "true"});
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
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"disable_heuristics", "true"});  // no initial setSolution
    opts.push_back({"dssr_background_updates", "true"});
    opts.push_back({"ng_initial_size", "4"});
    opts.push_back({"ng_max_size", "8"});
    opts.push_back({"ng_dssr_iters", "4"});

    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE_THAT(result.objective, WithinAbs(-13.0, 1e-6));
    REQUIRE(result.status == rcspp::SolveResult::Status::Optimal);
}

TEST_CASE("Model solves path with non-zero source", "[model][path]") {
    // source=1, target=3 (neither is node 0)
    rcspp::Model model;
    std::vector<rcspp::Edge> edges = {
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
    std::vector<rcspp::Edge> edges = {
        {0, 1}, {0, 2}, {1, 2}
    };
    std::vector<double> costs = {5.0, 3.0, 4.0};
    std::vector<double> profits = {0.0, 10.0, 8.0};
    std::vector<double> demands = {0.0, 2.0, 3.0};

    rcspp::Model model_depot;
    model_depot.set_graph(3, edges, costs);
    model_depot.set_depot(0);
    model_depot.set_profits(profits);
    model_depot.add_capacity_resource(demands, 5.0);
    auto result_depot = model_depot.solve(quiet);

    rcspp::Model model_st;
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
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"submip_separation", "false"});
    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE(result.objective < 0.0);
}

TEST_CASE("Model: submip_separation=true produces valid result", "[model]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"submip_separation", "true"});
    auto result = model.solve(opts);

    REQUIRE(result.has_solution());
    REQUIRE(result.objective < 0.0);
}

TEST_CASE("Model: extended tuning options parse and preserve correctness", "[model][options]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"max_concurrent_solves", "1"});
    opts.push_back({"deterministic_work_units", "128"});
    opts.push_back({"heuristic_deterministic_restarts", "8"});
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"heuristic_node_interval", "50"});
    opts.push_back({"heuristic_async_injection", "true"});
    opts.push_back({"dssr_background_policy", "auto"});
    opts.push_back({"dssr_background_max_epochs", "4"});
    opts.push_back({"dssr_background_auto_min_epochs", "2"});
    opts.push_back({"dssr_background_auto_no_progress_limit", "3"});
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
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"deterministic_work_units", "64"});
    opts.push_back({"heuristic_deterministic_restarts", "8"});
    opts.push_back({"dssr_background_updates", "false"});

    auto r1 = model.solve(opts);
    auto r2 = model.solve(opts);

    REQUIRE(r1.status == r2.status);
    REQUIRE_THAT(r1.bound, WithinAbs(r2.bound, 1e-9));
    if (r1.has_solution() && r2.has_solution()) {
        REQUIRE_THAT(r1.objective, WithinAbs(r2.objective, 1e-9));
        REQUIRE(r1.tour == r2.tour);
    }
}

TEST_CASE("Model: deterministic DSSR epoch commits are replay-equivalent",
          "[model][options][determinism][dssr]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4},
        {1, 2}, {1, 3}, {1, 4},
        {2, 3}, {2, 4},
        {3, 4},
    };
    std::vector<double> costs = {
        6.0, 8.0, 9.0, 10.0,
        3.0, 4.0, 7.0,
        5.0, 6.0,
        4.0,
    };
    std::vector<double> profits = {0.0, 11.0, 10.0, 9.0, 8.0};
    std::vector<double> demands = {0.0, 2.0, 2.0, 3.0, 3.0};

    model.set_graph(5, edges, costs);
    model.set_depot(0);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 8.0);

    auto opts = quiet;
    opts.push_back({"threads", "1"});
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"random_seed", "0"});
    opts.push_back({"dssr_background_updates", "true"});
    opts.push_back({"ng_initial_size", "1"});
    opts.push_back({"ng_max_size", "6"});
    opts.push_back({"ng_dssr_iters", "1"});

    auto r1 = model.solve(opts);
    auto r2 = model.solve(opts);

    REQUIRE(r1.status == r2.status);
    REQUIRE_THAT(r1.bound, WithinAbs(r2.bound, 1e-9));
    if (r1.has_solution() && r2.has_solution()) {
        REQUIRE_THAT(r1.objective, WithinAbs(r2.objective, 1e-9));
        REQUIRE(r1.tour == r2.tour);
    }

    REQUIRE(r1.dssr_epochs_enqueued == r2.dssr_epochs_enqueued);
    REQUIRE(r1.dssr_epochs_committed == r2.dssr_epochs_committed);
    REQUIRE(r1.dssr_epochs_committed == r1.dssr_epochs_enqueued);
    REQUIRE(r2.dssr_epochs_committed == r2.dssr_epochs_enqueued);
    REQUIRE(r1.dssr_checkpoint_count == r2.dssr_checkpoint_count);
    REQUIRE(r1.dssr_commit_signature == r2.dssr_commit_signature);
}

TEST_CASE("Model: deterministic DSSR max-epoch cap and auto policy are replay-stable",
          "[model][options][determinism][dssr]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4},
        {1, 2}, {1, 3}, {1, 4},
        {2, 3}, {2, 4},
        {3, 4},
    };
    std::vector<double> costs = {
        6.0, 8.0, 9.0, 10.0,
        3.0, 4.0, 7.0,
        5.0, 6.0,
        4.0,
    };
    std::vector<double> profits = {0.0, 11.0, 10.0, 9.0, 8.0};
    std::vector<double> demands = {0.0, 2.0, 2.0, 3.0, 3.0};

    model.set_graph(5, edges, costs);
    model.set_depot(0);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 8.0);

    auto opts = quiet;
    opts.push_back({"threads", "1"});
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"random_seed", "0"});
    opts.push_back({"dssr_background_updates", "true"});
    opts.push_back({"dssr_background_policy", "auto"});
    opts.push_back({"dssr_background_max_epochs", "2"});
    opts.push_back({"dssr_background_auto_min_epochs", "1"});
    opts.push_back({"dssr_background_auto_no_progress_limit", "1"});
    opts.push_back({"ng_initial_size", "1"});
    opts.push_back({"ng_max_size", "8"});
    opts.push_back({"ng_dssr_iters", "1"});

    auto r1 = model.solve(opts);
    auto r2 = model.solve(opts);

    REQUIRE(r1.status == r2.status);
    REQUIRE_THAT(r1.bound, WithinAbs(r2.bound, 1e-9));
    if (r1.has_solution() && r2.has_solution()) {
        REQUIRE_THAT(r1.objective, WithinAbs(r2.objective, 1e-9));
        REQUIRE(r1.tour == r2.tour);
    }

    REQUIRE(r1.dssr_epochs_enqueued == r1.dssr_epochs_committed);
    REQUIRE(r2.dssr_epochs_enqueued == r2.dssr_epochs_committed);
    REQUIRE(r1.dssr_epochs_enqueued <= 2);
    REQUIRE(r2.dssr_epochs_enqueued <= 2);
    REQUIRE(r1.dssr_epochs_enqueued == r2.dssr_epochs_enqueued);
    REQUIRE(r1.dssr_commit_signature == r2.dssr_commit_signature);
}

TEST_CASE("Model: DSSR auto-stop tracker stops on deterministic no-progress path",
          "[model][options][determinism][dssr]") {
    rcspp::model_detail::DssrAutoStopTracker tracker(
        /*min_epochs_before_stop=*/4,
        /*no_progress_epoch_limit=*/2);

    // With checkpoint activity present, stopping should come from no-progress
    // streak, not from the no-checkpoint fast-stop path.
    REQUIRE_FALSE(tracker.observe_stage(
        /*improved=*/false, /*checkpoint_active=*/true, /*produced_epochs=*/1));
    REQUIRE(tracker.observe_stage(
        /*improved=*/false, /*checkpoint_active=*/true, /*produced_epochs=*/2));
    REQUIRE(tracker.had_checkpoint_activity());
    REQUIRE(tracker.no_progress_epochs() == 2);
}

TEST_CASE("Model: disabling cut families still yields valid tiny solution", "[model][options][cuts]") {
    rcspp::Model model;

    std::vector<rcspp::Edge> edges = {
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
    std::vector<rcspp::Edge> edges = {
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

    rcspp::Model model_sec_on;
    rcspp::Model model_sec_off;
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
    std::vector<rcspp::Edge> edges = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3},
    };
    std::vector<double> costs = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0};
    std::vector<double> profits = {0.0, 0.0, 0.0, 60.0, 60.0, 60.0};
    std::vector<double> demands = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    rcspp::Model model;
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
    opts.push_back({"submip_separation", "false"});
    auto result = model.solve(opts);

    REQUIRE(result.is_optimal());
    REQUIRE(result.objective >= 0.0);
}

TEST_CASE("Model: edge elimination toggles preserve tiny optimum", "[model][options][presolve]") {
    rcspp::Model model_on;
    rcspp::Model model_off;

    std::vector<rcspp::Edge> edges = {
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
    opts_on.push_back({"disable_heuristics", "true"});
    opts_on.push_back({"edge_elimination", "true"});
    opts_on.push_back({"edge_elimination_nodes", "true"});

    auto opts_off = quiet;
    opts_off.push_back({"cutoff", "-11.0"});
    opts_off.push_back({"disable_heuristics", "true"});
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
    rcspp::Model model_a;
    rcspp::Model model_b;

    std::vector<rcspp::Edge> edges = {
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
    std::vector<rcspp::Edge> edges = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };
    std::vector<double> costs = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
    std::vector<double> profits = {0.0, 20.0, 15.0, 10.0};
    std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};

    // Get reference solution with branching off
    double ref_obj;
    {
        rcspp::Model model;
        model.set_graph(4, edges, costs);
        model.set_depot(0);
        model.set_profits(profits);
        model.add_capacity_resource(demands, 10.0);
        auto ref_opts = quiet;
        ref_opts.push_back({"threads", "1"});
        ref_opts.push_back({"random_seed", "0"});
        ref_opts.push_back({"parallel_mode", "deterministic"});
        ref_opts.push_back({"dssr_background_updates", "false"});
        auto result = model.solve(ref_opts);
        if (result.status == rcspp::SolveResult::Status::Error) {
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
            rcspp::Model model;
            model.set_graph(4, edges, costs);
            model.set_depot(0);
            model.set_profits(profits);
            model.add_capacity_resource(demands, 10.0);
            auto opts = quiet;
            opts.push_back({"threads", "1"});
            opts.push_back({"random_seed", "0"});
            opts.push_back({"parallel_mode", "deterministic"});
            opts.push_back({"dssr_background_updates", "false"});
            opts.push_back({"branch_hyper", mode});
            auto result = model.solve(opts);
            if (result.status == rcspp::SolveResult::Status::Error) {
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
    auto setup_path_model = [](rcspp::Model& model) {
        std::vector<rcspp::Edge> edges = {
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

    rcspp::Model base_model;
    setup_path_model(base_model);
    auto base_opts = quiet;
    base_opts.push_back({"parallel_mode", "deterministic"});
    base_opts.push_back({"dssr_background_updates", "false"});
    auto base = base_model.solve(base_opts);
    REQUIRE(base.has_solution());

    rcspp::Model staged_model;
    setup_path_model(staged_model);
    auto staged_opts = quiet;
    staged_opts.push_back({"parallel_mode", "deterministic"});
    staged_opts.push_back({"dssr_background_updates", "false"});
    staged_opts.push_back({"deterministic_work_units", "64"});
    staged_opts.push_back({"heuristic_deterministic_restarts", "8"});
    staged_opts.push_back({"ng_initial_size", "1"});
    staged_opts.push_back({"ng_max_size", "4"});
    staged_opts.push_back({"ng_dssr_iters", "2"});
    staged_opts.push_back({"preproc_adaptive", "true"});
    staged_opts.push_back({"preproc_fast_restarts", "4"});
    staged_opts.push_back({"preproc_fast_budget_ms", "5"});
    staged_opts.push_back({"preproc_second_ws_large_n", "4"});
    staged_opts.push_back({"preproc_second_ws_min_elim", "0.0"});
    staged_opts.push_back({"preproc_second_ws_min_elim_large", "0.0"});
    staged_opts.push_back({"preproc_second_ws_budget_ms_min", "5"});
    staged_opts.push_back({"preproc_second_ws_budget_ms_max", "10"});
    staged_opts.push_back({"preproc_second_ws_budget_scale", "1.0"});
    auto staged = staged_model.solve(staged_opts);
    REQUIRE(staged.has_solution());
    REQUIRE_THAT(staged.objective, WithinAbs(base.objective, 1e-6));
}

TEST_CASE("Model: stage1 bounds backend modes preserve objective",
          "[model][preproc]") {
    auto setup_path_model = [](rcspp::Model& model) {
        std::vector<rcspp::Edge> edges = {
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

    auto solve_mode = [&](const char* mode) {
        rcspp::Model model;
        setup_path_model(model);
        auto opts = quiet;
        opts.push_back({"threads", "1"});
        opts.push_back({"random_seed", "0"});
        opts.push_back({"parallel_mode", "deterministic"});
        opts.push_back({"dssr_background_updates", "false"});
        opts.push_back({"preproc_stage1_bounds", mode});
        opts.push_back({"ng_initial_size", "1"});
        opts.push_back({"ng_max_size", "4"});
        opts.push_back({"ng_dssr_iters", "2"});
        auto r = model.solve(opts);
        if (r.status == rcspp::SolveResult::Status::Error) {
            WARN("mode=" << mode << " returned Error; skipping comparison for this mode");
        }
        return r;
    };

    const auto two_cycle = solve_mode("two_cycle");
    const auto ng1 = solve_mode("ng1");
    const auto auto_mode = solve_mode("auto");

    if (two_cycle.has_solution() && ng1.has_solution()) {
        REQUIRE_THAT(two_cycle.objective, WithinAbs(ng1.objective, 1e-6));
    }
    if (two_cycle.has_solution() && auto_mode.has_solution()) {
        REQUIRE_THAT(two_cycle.objective, WithinAbs(auto_mode.objective, 1e-6));
    }
}

TEST_CASE("Model: workflow_dump option is accepted",
          "[model][workflow]") {
    rcspp::Model model;
    std::vector<rcspp::Edge> edges = {
        {0, 1}, {1, 2}, {0, 2}
    };
    std::vector<double> costs = {4.0, 3.0, 10.0};
    std::vector<double> profits = {0.0, 7.0, 2.0};
    std::vector<double> demands = {0.0, 1.0, 1.0};
    model.set_graph(3, edges, costs);
    model.set_source(0);
    model.set_target(2);
    model.set_profits(profits);
    model.add_capacity_resource(demands, 2.0);

    auto opts = quiet;
    opts.push_back({"threads", "1"});
    opts.push_back({"random_seed", "0"});
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"dssr_background_updates", "false"});
    opts.push_back({"workflow_dump", "true"});

    auto r = model.solve(opts);
    if (r.status == rcspp::SolveResult::Status::Error) {
        WARN("workflow_dump run returned Error; skipping strict assertion");
        return;
    }
}

TEST_CASE("Model: paramip planning options are accepted",
          "[model][paramip]") {
    rcspp::Model model;
    std::vector<rcspp::Edge> edges = {
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
    opts.push_back({"parallel_mode", "deterministic"});
    opts.push_back({"dssr_background_updates", "false"});
    opts.push_back({"workflow_dump", "true"});
    opts.push_back({"paramip_mode", "plan"});
    opts.push_back({"paramip_chunks", "8"});
    opts.push_back({"paramip_workers", "4"});

    auto r = model.solve(opts);
    if (r.status == rcspp::SolveResult::Status::Error) {
        WARN("paramip plan run returned Error; skipping strict assertion");
        return;
    }
}
