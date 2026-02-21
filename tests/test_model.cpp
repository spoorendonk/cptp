#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "model/model.h"

using Catch::Matchers::WithinAbs;

static const cptp::SolverOptions quiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

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
