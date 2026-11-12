// m capacity rows end to end (#90): build, loaders, heuristic residuals, and
// solves where a solution must satisfy every row.
#include "core/io.h"
#include "core/problem.h"
#include "heuristic/primal_heuristic.h"
#include "model/model.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;

namespace {

const cptp::SolverOptions kQuiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

// Complete graph on depot 0 and customers 1..3, every edge costs 1, so a tour
// through j customers costs j + 1 (one customer: depot edge twice, cost 2).
// Profits 12, 10, 11. The rows:
//   row 0: d = (0, 1, 1, 1), Q = 3 -> all three fit: 33 - 4 = -29
//   row 1: d = (0, 2, 2, 1), Q = 3 -> {1, 3} or {2, 3}: {1, 3} gives 23 - 3 = -20
//   row 2: d = (0, 1, 0, 1), Q = 1 -> not both 1 and 3: {2, 3} gives 21 - 3 = -18
const std::vector<cptp::Edge> kEdges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
const std::vector<double> kCosts(6, 1.0);
const std::vector<double> kProfits = {0.0, 12.0, 10.0, 11.0};
const std::vector<std::vector<double>> kRows = {
    {0.0, 1.0, 1.0, 1.0}, {0.0, 2.0, 2.0, 1.0}, {0.0, 1.0, 0.0, 1.0}};
const std::vector<double> kCaps = {3.0, 3.0, 1.0};

void fill_model(cptp::Model& model, int32_t rows) {
    model.set_graph(4, kEdges, kCosts);
    model.set_profits(kProfits);
    for (int32_t k = 0; k < rows; ++k) {
        model.add_capacity_resource(kRows[k], kCaps[k]);
    }
}

cptp::Problem make_problem(int32_t rows) {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits,
               std::vector<std::vector<double>>(kRows.begin(), kRows.begin() + rows),
               std::vector<double>(kCaps.begin(), kCaps.begin() + rows), 0, 0);
    return prob;
}

// Every row's load of the visited nodes is within its capacity.
void require_tour_fits(const cptp::Problem& prob, const std::vector<int32_t>& tour) {
    std::vector<bool> seen(static_cast<size_t>(prob.num_nodes()), false);
    for (int32_t v : tour) {
        seen[static_cast<size_t>(v)] = true;
    }
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        double load = 0.0;
        for (int32_t v = 0; v < prob.num_nodes(); ++v) {
            if (seen[static_cast<size_t>(v)]) {
                load += prob.demand(k, v);
            }
        }
        INFO("row " << k);
        REQUIRE(load <= prob.capacity(k) + 1e-9);
    }
}

std::filesystem::path write_temp(const std::string& name, const std::string& text) {
    auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream(path) << text;
    return path;
}

}  // namespace

TEST_CASE("Multicap: Problem keeps every row and row 0 as the scalar view", "[multicap]") {
    const auto prob = make_problem(3);
    REQUIRE(prob.num_capacities() == 3);
    REQUIRE(prob.capacities() == kCaps);
    REQUIRE(prob.demand_rows() == kRows);
    REQUIRE(prob.capacity() == kCaps[0]);
    REQUIRE(prob.demand(2) == kRows[0][2]);
    REQUIRE(prob.demand(1, 2) == kRows[1][2]);
    REQUIRE(prob.demands(2) == kRows[2]);
    REQUIRE(prob.fits_alone(1));

    // A customer whose demand exceeds one row's capacity does not fit alone.
    cptp::Problem tight;
    tight.build(4, kEdges, kCosts, kProfits, {{0.0, 1.0, 1.0, 1.0}, {0.0, 1.0, 5.0, 1.0}},
                {3.0, 3.0}, 0, 0);
    REQUIRE(tight.fits_alone(1));
    REQUIRE_FALSE(tight.fits_alone(2));
}

TEST_CASE("Multicap: the optimum satisfies every row", "[multicap][model]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        INFO("m = " << c.rows);
        cptp::Model model;
        fill_model(model, c.rows);
        model.set_depot(0);
        const auto result = model.solve(kQuiet);
        REQUIRE(result.is_optimal());
        REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
        require_tour_fits(model.problem(), result.tour);
    }
}

TEST_CASE("Multicap: Q_k = 0 is a limit, Q_k < 0 is rejected", "[multicap][model]") {
    // Row 1 has Q = 0 and only customer 1 has demand on it, so 1 is out and
    // the optimum is {2, 3}: 21 - 3 = -18 (row 0 alone would allow all three).
    cptp::Model model;
    model.set_graph(4, kEdges, kCosts);
    model.set_profits(kProfits);
    model.set_depot(0);
    model.add_capacity_resource(kRows[0], kCaps[0]);
    model.add_capacity_resource(std::vector<double>{0.0, 1.0, 0.0, 0.0}, 0.0);
    const auto result = model.solve(kQuiet);
    REQUIRE(result.is_optimal());
    REQUIRE_THAT(result.objective, WithinAbs(-18.0, 1e-6));
    require_tour_fits(model.problem(), result.tour);

    cptp::Problem bad;
    REQUIRE_THROWS_AS(
        bad.build(4, kEdges, kCosts, kProfits, {kRows[0], kRows[1]}, {3.0, -1.0}, 0, 0),
        std::invalid_argument);
}

TEST_CASE("Multicap: s-t path respects every row", "[multicap][model][path]") {
    // Path 0 -> 3: endpoints count against every row. Row 1 leaves room for
    // one of customers 1 and 2 (d = 2 each, Q - d_3 = 2); row 2 forbids 1.
    cptp::Model model;
    fill_model(model, 3);
    model.set_source(0);
    model.set_target(3);
    const auto result = model.solve(kQuiet);
    REQUIRE(result.is_optimal());
    // 0 -> 2 -> 3: cost 2, profit 10 + 11.
    REQUIRE_THAT(result.objective, WithinAbs(-19.0, 1e-6));
    require_tour_fits(model.problem(), result.tour);
}

TEST_CASE("Multicap: Residual tests and updates every row", "[multicap][heuristic]") {
    const auto prob = make_problem(3);
    cptp::heuristic::detail::Residual rem(prob);
    REQUIRE(rem.size() == 3);
    REQUIRE(rem[0] == 3.0);
    REQUIRE(rem[2] == 1.0);

    REQUIRE(rem.fits(prob, 1));
    rem.take(prob, 1);  // (2, 1, 0)
    REQUIRE(rem[1] == 1.0);
    REQUIRE(rem[2] == 0.0);
    REQUIRE_FALSE(rem.fits(prob, 3));          // row 2: d = 1 > 0
    REQUIRE_FALSE(rem.fits(prob, 2));          // row 1: d = 2 > 1
    REQUIRE(rem.fits_swap(prob, 1, 3, 1e-9));  // (2, 2, 0)
    REQUIRE(rem.at_least(0.0));
    REQUIRE_FALSE(rem.at_least(0.5));

    rem.swap(prob, 1, 3);
    REQUIRE(rem[0] == 2.0);
    REQUIRE(rem[1] == 2.0);
    REQUIRE(rem[2] == 0.0);
    rem.give(prob, 3);
    REQUIRE(rem[0] == 3.0);
    REQUIRE(rem[1] == 3.0);
    REQUIRE(rem[2] == 1.0);

    // Tour {2, 3} leaves (1, 0, 0): swapping 2 for 1 fits rows 0 and 1 but
    // not row 2 (0 + 0 - 1 < 0).
    cptp::heuristic::detail::Residual pair(prob);
    pair.take(prob, 2);
    pair.take(prob, 3);
    REQUIRE_FALSE(pair.fits_swap(prob, 2, 1, 1e-9));
}

TEST_CASE("Multicap: initial heuristic solution satisfies every row", "[multicap][heuristic]") {
    const auto prob = make_problem(3);
    const auto start = cptp::heuristic::build_initial_solution(prob, 8);
    REQUIRE(start.objective < std::numeric_limits<double>::max());
    const int32_t m = prob.num_edges();
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        double load = 0.0;
        for (int32_t v = 0; v < prob.num_nodes(); ++v) {
            load += prob.demand(k, v) *
                    start.col_values[static_cast<size_t>(m) + static_cast<size_t>(v)];
        }
        INFO("row " << k);
        REQUIRE(load <= prob.capacity(k) + 1e-9);
    }
    REQUIRE(start.objective >= -18.0 - 1e-9);
}

TEST_CASE("Multicap: TSPLIB loader reads m capacities and m demand columns", "[multicap][io]") {
    const auto prob =
        cptp::io::load("benchmarks/instances/multicap/spprclib/m2/B-n45-k6-54-m2.sppcc");
    REQUIRE(prob.num_capacities() == 2);
    REQUIRE(prob.capacity(0) == 100.0);
    REQUIRE(prob.capacity(1) == 112.0);
    REQUIRE(prob.num_nodes() == 45);
    REQUIRE(prob.demand(0, prob.source()) == 0.0);
    REQUIRE(prob.demand(1, prob.source()) == 0.0);

    const auto bad = write_temp("cptp_multicap_bad.vrp",
                                "NAME : bad\nTYPE : CVRP\nDIMENSION : 3\n"
                                "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 5 7 9\n"
                                "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                "DEMAND_SECTION\n1 0 0\n2 2 3\n3 1 4\n"
                                "DEPOT_SECTION\n1\n-1\nEOF\n");
    REQUIRE_THROWS_AS(cptp::io::load_tsplib(bad), std::runtime_error);
    std::filesystem::remove(bad);

    // No CAPACITY: one row with no limit, not a limit of 0.
    const auto open = write_temp("cptp_multicap_open.vrp",
                                 "NAME : open\nTYPE : CVRP\nDIMENSION : 3\n"
                                 "EDGE_WEIGHT_TYPE : EUC_2D\n"
                                 "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                 "DEPOT_SECTION\n1\n-1\nEOF\n");
    const auto unlimited = cptp::io::load_tsplib(open);
    std::filesystem::remove(open);
    REQUIRE(unlimited.num_capacities() == 1);
    REQUIRE(unlimited.capacity() >= 1e17);
}

TEST_CASE("Multicap: numeric loader reads m rows", "[multicap][io]") {
    // 3 nodes, 6 arcs (one kept per edge), m = 2, then capacities, then s t.
    const auto path = write_temp("cptp_multicap.txt",
                                 "3 6 2\n"
                                 "0 0 0 0\n"
                                 "1 5 2 3\n"
                                 "2 4 1 6\n"
                                 "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n"
                                 "4 7\n"
                                 "0 0\n");
    const auto prob = cptp::io::load_numeric(path);
    std::filesystem::remove(path);
    REQUIRE(prob.num_capacities() == 2);
    REQUIRE(prob.capacities() == std::vector<double>{4.0, 7.0});
    REQUIRE(prob.demands(1) == std::vector<double>{0.0, 3.0, 6.0});
    REQUIRE(prob.num_edges() == 3);

    const auto short_caps =
        write_temp("cptp_multicap_short.txt", "2 2 2\n0 0 0 0\n1 5 2 3\n0 1 1\n1 0 1\n4\n");
    REQUIRE_THROWS_AS(cptp::io::load_numeric(short_caps), std::runtime_error);
    std::filesystem::remove(short_caps);
}
