// m capacity rows end to end (#90): build, loaders, heuristic residuals, and
// solves where a solution must satisfy every row.
#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/io.h"
#include "core/problem.h"
#include "heuristic/primal_heuristic.h"
#include "model/hyperplane.h"
#include "model/model.h"
#include "preprocess/edge_elimination.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/rglm_separator.h"
#include "sep/separation_context.h"
#include "sep/separator.h"

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
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

    // Trailing text after the capacity is ignored, as it was for m = 1.
    const auto trailing = write_temp("cptp_multicap_trailing.vrp",
                                     "NAME : trailing\nTYPE : CVRP\nDIMENSION : 3\n"
                                     "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 100 (vehicles)\n"
                                     "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                     "DEMAND_SECTION\n1 0\n2 2\n3 1\n"
                                     "DEPOT_SECTION\n1\n-1\nEOF\n");
    const auto one_row = cptp::io::load_tsplib(trailing);
    std::filesystem::remove(trailing);
    REQUIRE(one_row.num_capacities() == 1);
    REQUIRE(one_row.capacity() == 100.0);

    // A CAPACITY without a leading number is still an error.
    const auto none = write_temp("cptp_multicap_none.vrp",
                                 "NAME : none\nTYPE : CVRP\nDIMENSION : 3\n"
                                 "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : many\n"
                                 "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                 "DEPOT_SECTION\n1\n-1\nEOF\n");
    REQUIRE_THROWS_AS(cptp::io::load_tsplib(none), std::runtime_error);
    std::filesystem::remove(none);
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

    // With m > 1 a missing capacity line is an error, not "no limit".
    const auto no_caps =
        write_temp("cptp_multicap_nocaps.txt", "2 2 2\n0 0 0 0\n1 5 2 3\n0 1 1\n1 0 1\n");
    REQUIRE_THROWS_AS(cptp::io::load_numeric(no_caps), std::runtime_error);
    std::filesystem::remove(no_caps);

    // Truncated node lines are an error.
    const auto truncated = write_temp("cptp_multicap_truncated.txt", "3 0 2\n0 0 0 0\n1 5 2\n");
    REQUIRE_THROWS_AS(cptp::io::load_numeric(truncated), std::runtime_error);
    std::filesystem::remove(truncated);
}

namespace {

// A point where every customer is served and x(delta({1, 2, 3})) = 1 is the
// smallest cut, so {1, 2, 3} is a Gomory-Hu candidate: x01 = x03 = x13 = 0.5,
// x12 = x23 = 1 (each customer has degree 2). make_point(prob, t) scales x and
// the customers' y by t, so each customer keeps degree 2 y_i.
struct LpPoint {
    std::vector<double> x;
    std::vector<double> y{1.0, 1.0, 1.0, 1.0};
    std::unique_ptr<cptp::gomory_hu_tree> tree;
};

double point_value(int32_t u, int32_t v) {
    if ((u == 1 && v == 2) || (u == 2 && v == 3)) {
        return 1.0;
    }
    if ((u == 0 && v == 1) || (u == 0 && v == 3) || (u == 1 && v == 3)) {
        return 0.5;
    }
    return 0.0;
}

LpPoint make_point(const cptp::Problem& prob, double scale = 1.0) {
    LpPoint p;
    for (size_t i = 1; i < p.y.size(); ++i) {
        p.y[i] = scale;
    }
    const auto& g = prob.graph();
    p.x.assign(static_cast<size_t>(prob.num_edges()), 0.0);
    cptp::digraph_builder builder(prob.num_nodes());
    for (auto e : g.edges()) {
        const int32_t u = g.edge_source(e);
        const int32_t v = g.edge_target(e);
        const double xe = scale * point_value(std::min(u, v), std::max(u, v));
        if (xe > 0.0) {
            p.x[static_cast<size_t>(e)] = xe;
            builder.add_arc(u, v, xe);
            builder.add_arc(v, u, xe);
        }
    }
    auto [sg, cap] = builder.build();
    p.tree = std::make_unique<cptp::gomory_hu_tree>(sg, cap, prob.source());
    return p;
}

std::vector<cptp::sep::Cut> cuts_of(const cptp::Problem& prob, cptp::sep::Separator& separator,
                                    double scale = 1.0) {
    const auto point = make_point(prob, scale);
    const cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = point.x,
        .y_values = point.y,
        .x_offset = 0,
        .y_offset = prob.num_edges(),
        .tol = 1e-6,
        .flow_tree = point.tree.get(),
    };
    return separator.separate(ctx);
}

std::vector<cptp::sep::Cut> rci_cuts(const cptp::Problem& prob, cptp::sep::CapacityRows rows) {
    cptp::sep::RCISeparator rci(rows);
    return cuts_of(prob, rci);
}

// Row 0 loose (d = 1 each, Q = 10), row 1 tight (d = 2, 2, 1, Q = 3).
cptp::Problem make_loose_first_row() {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, std::vector<std::vector<double>>{kRows[0], kRows[1]},
               std::vector<double>{10.0, kCaps[1]}, 0, 0);
    return prob;
}

// Row 1 of the instance as the only row.
cptp::Problem make_single_tight_row() {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, std::vector<std::vector<double>>{kRows[1]},
               std::vector<double>{kCaps[1]}, 0, 0);
    return prob;
}

void require_same_cuts(const std::vector<cptp::sep::Cut>& a, const std::vector<cptp::sep::Cut>& b) {
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].indices == b[i].indices);
        REQUIRE(a[i].values == b[i].values);
        REQUIRE(a[i].rhs == b[i].rhs);
    }
}

}  // namespace

TEST_CASE("Multicap: RCI separates from rows after the first", "[multicap][rci]") {
    using cptp::sep::CapacityRows;
    // Row 0 (d = 1 each, Q = 3) gives no RCI: no S has d(S) > Q. Row 1
    // (d = 2, 2, 1, Q = 3) cuts S = {1, 2, 3}: d(S) = 5 needs two vehicles.
    const auto one_row = make_problem(1);
    REQUIRE(rci_cuts(one_row, CapacityRows::kAll).empty());

    const auto two_rows = make_problem(2);
    REQUIRE(rci_cuts(two_rows, CapacityRows::kFirst).empty());
    const auto all = rci_cuts(two_rows, CapacityRows::kAll);
    const auto best = rci_cuts(two_rows, CapacityRows::kMostViolated);
    REQUIRE_FALSE(all.empty());
    REQUIRE_FALSE(best.empty());
    REQUIRE(best.size() <= all.size());
    for (const auto& cut : all) {
        REQUIRE(cut.violation > 1e-6);
    }
}

TEST_CASE("Multicap: with one row RCI is the same for every row option", "[multicap][rci]") {
    using cptp::sep::CapacityRows;
    // Row 1 of the instance as the only row, so RCI has cuts to give.
    const auto prob = make_single_tight_row();
    const auto all = rci_cuts(prob, CapacityRows::kAll);
    REQUIRE_FALSE(all.empty());
    for (const auto rows : {CapacityRows::kMostViolated, CapacityRows::kFirst}) {
        require_same_cuts(rci_cuts(prob, rows), all);
    }
}

TEST_CASE("Multicap: capacity_cut_rows takes all, most_violated or first", "[multicap][rci]") {
    for (const char* value : {"all", "most_violated", "first"}) {
        cptp::Model model;
        fill_model(model, 3);
        auto options = kQuiet;
        options.emplace_back("capacity_cut_rows", value);
        const auto result = model.solve(options);
        INFO(value);
        REQUIRE(result.is_optimal());
        REQUIRE_THAT(result.objective, WithinAbs(-18.0, 1e-6));
    }
    REQUIRE_THROWS_AS(cptp::sep::parse_capacity_rows("some"), std::invalid_argument);
}

TEMPLATE_TEST_CASE("Multicap: GLM-type cuts come from rows after the first", "[multicap][glm]",
                   cptp::sep::MultistarSeparator, cptp::sep::RGLMSeparator) {
    using cptp::sep::CapacityRows;
    // At the point, S = {1, 2, 3} with x(delta(S)) = 1. Row 0 with Q = 10
    // gives no violated multistar (0.6 < 1) or RGLM (ceil(3/10) = 1); row 1
    // (d(S) = 5, Q = 3) gives both.
    const auto prob = make_loose_first_row();
    TestType first(CapacityRows::kFirst);
    TestType all(CapacityRows::kAll);
    TestType best(CapacityRows::kMostViolated);
    REQUIRE(cuts_of(prob, first).empty());
    const auto all_cuts = cuts_of(prob, all);
    REQUIRE_FALSE(all_cuts.empty());
    require_same_cuts(cuts_of(prob, best), all_cuts);  // one row violated: same cuts

    // With one row, every option gives the same cuts.
    const auto single = make_single_tight_row();
    TestType single_all(CapacityRows::kAll);
    const auto reference = cuts_of(single, single_all);
    REQUIRE_FALSE(reference.empty());
    for (const auto rows : {CapacityRows::kMostViolated, CapacityRows::kFirst}) {
        TestType other(rows);
        require_same_cuts(cuts_of(single, other), reference);
    }
}

namespace {

// The example of the #86 write-up (#105) on S = {1, 2, 3}: row 0 has
// d = 7 each and Q = 10 (rho = 2.1), row 1 has d = 3 each and Q = 10
// (rho = 0.9). At the point scaled by t = 10/21 (y_i = t, so row 0's served
// load is exactly 10 and both rows hold), row 0's cut has
// 2(3 - (1 - t) 2.1 / 0.1) < 0 on the right and row 1 has rho < 1, so
// neither row gives a violated RCI or RGLM. The surrogate row has s_i = 1
// each against 2: rho = 1.5, remainder 1, rhs 2(2 - 3(1 - t)) = 6t - 2, and
// x(delta(S)) = t, so it is violated by 6t - 2 - t = 8/21.
constexpr double kSurrogateScale = 10.0 / 21.0;

cptp::Problem make_surrogate_example() {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits,
               std::vector<std::vector<double>>{{0.0, 7.0, 7.0, 7.0}, {0.0, 3.0, 3.0, 3.0}},
               std::vector<double>{10.0, 10.0}, 0, 0);
    return prob;
}

}  // namespace

TEMPLATE_TEST_CASE("Multicap: the surrogate row cuts where no single row does", "[multicap][rci]",
                   cptp::sep::RCISeparator, cptp::sep::RGLMSeparator) {
    using cptp::sep::CapacityRows;
    const auto prob = make_surrogate_example();
    for (const auto rows : {CapacityRows::kAll, CapacityRows::kMostViolated}) {
        TestType plain(rows);
        REQUIRE(cuts_of(prob, plain, kSurrogateScale).empty());
        TestType with_surrogate(rows, true);
        const auto cuts = cuts_of(prob, with_surrogate, kSurrogateScale);
        REQUIRE_FALSE(cuts.empty());
        double best = 0.0;
        for (const auto& cut : cuts) {
            best = std::max(best, cut.violation);
        }
        REQUIRE_THAT(best, WithinAbs(8.0 / 21.0, 1e-6));
    }

    // With one bounded row there is no surrogate: the same cuts as without.
    const auto single = make_single_tight_row();
    TestType plain(CapacityRows::kAll);
    TestType with_surrogate(CapacityRows::kAll, true);
    const auto reference = cuts_of(single, plain);
    REQUIRE_FALSE(reference.empty());
    require_same_cuts(cuts_of(single, with_surrogate), reference);
}

TEMPLATE_TEST_CASE("Multicap: under first the surrogate's cuts come next to row 0's",
                   "[multicap][rci]", cptp::sep::RCISeparator, cptp::sep::RGLMSeparator) {
    using cptp::sep::CapacityRows;
    // The example at t = 0.9: row 0 now cuts S too (by 41t - 36 for RCI) and
    // the surrogate by more (5t - 2). first keeps row 0's cuts and adds the
    // surrogate's; most_violated keeps one per candidate.
    constexpr double kScale = 0.9;
    const auto prob = make_surrogate_example();
    TestType first(CapacityRows::kFirst);
    const auto row0 = cuts_of(prob, first, kScale);
    REQUIRE_FALSE(row0.empty());
    TestType first_surrogate(CapacityRows::kFirst, true);
    const auto both = cuts_of(prob, first_surrogate, kScale);
    REQUIRE(both.size() > row0.size());
    for (const auto& cut : row0) {
        const bool kept = std::ranges::any_of(both, [&](const cptp::sep::Cut& other) {
            return other.indices == cut.indices && other.values == cut.values &&
                   other.rhs == cut.rhs;
        });
        REQUIRE(kept);
    }
    TestType best_surrogate(CapacityRows::kMostViolated, true);
    REQUIRE(cuts_of(prob, best_surrogate, kScale).size() < both.size());
}

TEST_CASE("Multicap: capacity_cut_surrogate keeps every row's optimum", "[multicap][rci]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        for (const char* rows : {"all", "most_violated", "first"}) {
            INFO("m = " << c.rows << ", capacity_cut_rows = " << rows);
            cptp::Model model;
            fill_model(model, c.rows);
            model.set_depot(0);
            auto options = kQuiet;
            options.emplace_back("capacity_cut_surrogate", "true");
            options.emplace_back("capacity_cut_rows", rows);
            options.emplace_back("enable_rglm", "true");
            const auto result = model.solve(options);
            REQUIRE(result.is_optimal());
            REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
            require_tour_fits(model.problem(), result.tour);
        }
    }
}

TEST_CASE("Multicap: multistar_largest_side keeps every row's optimum", "[multicap][glm]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        INFO("m = " << c.rows);
        cptp::Model model;
        fill_model(model, c.rows);
        model.set_depot(0);
        auto options = kQuiet;
        options.emplace_back("multistar_largest_side", "true");
        const auto result = model.solve(options);
        REQUIRE(result.is_optimal());
        REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
        require_tour_fits(model.problem(), result.tour);
    }
}

TEST_CASE("Multicap: enable_xcover keeps every row's optimum", "[multicap][xcover]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        INFO("m = " << c.rows);
        cptp::Model model;
        fill_model(model, c.rows);
        model.set_depot(0);
        auto options = kQuiet;
        options.emplace_back("enable_xcover", "true");
        const auto result = model.solve(options);
        REQUIRE(result.is_optimal());
        REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
        require_tour_fits(model.problem(), result.tour);
    }
    cptp::Model path;
    fill_model(path, 3);
    path.set_source(0);
    path.set_target(3);
    auto options = kQuiet;
    options.emplace_back("enable_xcover", "true");
    const auto result = path.solve(options);
    REQUIRE(result.is_optimal());
    REQUIRE_THAT(result.objective, WithinAbs(-19.0, 1e-6));
    require_tour_fits(path.problem(), result.tour);
}

TEST_CASE("Multicap: enable_xcover agrees with off on random two-row instances",
          "[multicap][xcover]") {
    // Random tours (t = 0) and s-t paths (t = 7) with two binding rows: the
    // cross-row cover cuts must not cut off the optimum found without them.
    std::mt19937 rng(1661);
    std::uniform_int_distribution<int> cost(1, 20);
    std::uniform_int_distribution<int> profit(0, 30);
    std::uniform_int_distribution<int> dem(0, 4);
    constexpr int32_t kNodes = 8;
    for (int trial = 0; trial < 16; ++trial) {
        std::vector<cptp::Edge> edges;
        std::vector<double> costs;
        for (int32_t u = 0; u < kNodes; ++u) {
            for (int32_t v = u + 1; v < kNodes; ++v) {
                edges.push_back({u, v});
                costs.push_back(cost(rng));
            }
        }
        std::vector<double> profits(kNodes);
        std::vector<std::vector<double>> rows(2, std::vector<double>(kNodes));
        for (int32_t v = 0; v < kNodes; ++v) {
            profits[static_cast<size_t>(v)] = profit(rng);
            for (auto& row : rows) {
                row[static_cast<size_t>(v)] = dem(rng);
            }
        }
        const int32_t target = trial % 2 == 0 ? 0 : kNodes - 1;
        cptp::Problem prob;
        prob.build(kNodes, edges, costs, profits, rows, std::vector<double>{8.0, 8.0}, 0, target);
        auto solve = [&](const char* xcover) {
            cptp::Model model;
            model.set_problem(prob);
            auto opts = kQuiet;
            opts.emplace_back("enable_xcover", xcover);
            const auto result = model.solve(opts);
            REQUIRE(result.is_optimal());
            require_tour_fits(model.problem(), result.tour);
            return result.objective;
        };
        INFO("trial " << trial);
        REQUIRE_THAT(solve("true"), WithinAbs(solve("false"), 1e-6));
    }
}

TEST_CASE("Multicap: most_violated keeps the most violated row's multistar", "[multicap][glm]") {
    using cptp::sep::CapacityRows;
    // Both rows of make_problem(2) cut S = {1, 2, 3}: row 0 by 2 - 1 = 1,
    // row 1 by 10/3 - 1 = 7/3. all keeps both per target, most_violated row 1's.
    const auto prob = make_problem(2);
    cptp::sep::MultistarSeparator all(CapacityRows::kAll);
    cptp::sep::MultistarSeparator best(CapacityRows::kMostViolated);
    const auto all_cuts = cuts_of(prob, all);
    const auto best_cuts = cuts_of(prob, best);
    REQUIRE_FALSE(best_cuts.empty());
    REQUIRE(all_cuts.size() == 2 * best_cuts.size());
    for (const auto& cut : best_cuts) {
        REQUIRE_THAT(cut.violation, WithinAbs(7.0 / 3.0, 1e-9));
    }
}

namespace {

// The hyperplane columns of `mode` (y columns from offset 6).
constexpr int32_t kYOffset = 6;

std::vector<cptp::detail::Hyperplane> hyperplanes_of(const cptp::Problem& prob,
                                                     const std::string& mode) {
    return cptp::detail::HyperplaneCandidates(prob, kYOffset, mode).columns();
}

std::vector<std::string> names_of(const std::vector<cptp::detail::Hyperplane>& hs) {
    std::vector<std::string> names;
    names.reserve(hs.size());
    for (const auto& h : hs) {
        names.push_back(h.name);
    }
    return names;
}

}  // namespace

TEST_CASE("Multicap: demand hyperplanes come from every bounded row", "[multicap][branching]") {
    const auto prob = make_problem(2);
    // One cluster {1, 2, 3}, weighed by each row's demands. Each family on its
    // own: under "all", row 0's unit demands make cluster_0, demand and
    // cardinality one hyperplane, kept once.
    REQUIRE(names_of(hyperplanes_of(prob, "clusters")) ==
            std::vector<std::string>{"cluster_0", "cluster_0_row1"});
    const auto hs = hyperplanes_of(prob, "demand");
    REQUIRE(names_of(hs) == std::vector<std::string>{"demand", "demand_row1"});
    // Under "all" demand_row1 is cluster_0_row1 again, and the one pair is (1, 2).
    REQUIRE(names_of(hyperplanes_of(prob, "all")) ==
            std::vector<std::string>{"pair_1_2", "cluster_0", "cluster_0_row1"});
    const auto row1 =
        std::ranges::find(hs, std::string("demand_row1"), &cptp::detail::Hyperplane::name);
    REQUIRE(row1 != hs.end());
    REQUIRE(row1->indices == std::vector<int32_t>{7, 8, 9});
    REQUIRE(row1->values == std::vector<double>{2.0, 2.0, 1.0});
}

TEST_CASE("Multicap: one row keeps the one-row hyperplane names", "[multicap][branching]") {
    const auto prob = make_problem(1);
    const auto demand = hyperplanes_of(prob, "demand");
    REQUIRE(names_of(demand) == std::vector<std::string>{"demand"});
    REQUIRE(demand[0].indices == std::vector<int32_t>{7, 8, 9});
    REQUIRE(demand[0].values == std::vector<double>{1.0, 1.0, 1.0});
    REQUIRE(names_of(hyperplanes_of(prob, "clusters")) == std::vector<std::string>{"cluster_0"});
    for (const auto& name : names_of(hyperplanes_of(prob, "all"))) {
        REQUIRE_FALSE(name.contains("_row"));
    }
}

TEST_CASE("Multicap: rows with Q_k = 0 or no limit give no hyperplane", "[multicap][branching]") {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, {kRows[0], kRows[1], kRows[2], kRows[1]},
               {3.0, 0.0, cptp::sep::kInfiniteCapacity, 3.0}, 0, 0);
    const cptp::detail::HyperplaneCandidates gen(prob, kYOffset, "demand");
    REQUIRE(gen.rows() == std::vector<int32_t>{0, 3});
    REQUIRE(names_of(hyperplanes_of(prob, "demand")) ==
            std::vector<std::string>{"demand", "demand_row3"});
}

TEST_CASE("Multicap: labeling over all rows fixes an edge row 0 alone does not",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    // Tour from depot 0 over customers 1 (profit 5) and 2 (profit 1); edges
    // (0,1) cost 10, (0,2) and (1,2) cost 1. Row 0 (d = 1, 1; Q = 2) admits
    // both, row 1 (d = 2, 2; Q = 3) only one, and each fits alone, so demand
    // reachability fixes nothing. The best tour is empty (0), so UB = 0.
    // Row 0 reaches 1 cheaply via 2: f[1] = 2 - 6 = -4, f[2] = 0, and edge
    // (1,2) has lb = -4 + 1 + 0 = -3. Row 1 must go direct: f[1] = 10 - 5 = 5,
    // lb = 5 + 1 + 0 = 6 > 0.
    cptp::Problem prob;
    prob.build(3, std::vector<cptp::Edge>{{0, 1}, {0, 2}, {1, 2}},
               std::vector<double>{10.0, 1.0, 1.0}, std::vector<double>{0.0, 5.0, 1.0},
               std::vector<std::vector<double>>{{0.0, 1.0, 1.0}, {0.0, 2.0, 2.0}},
               std::vector<double>{2.0, 3.0}, 0, 0);
    const auto first = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kFirst);
    const auto all = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kAll);
    REQUIRE_THAT(first[1], WithinAbs(-4.0, 1e-9));
    REQUIRE_THAT(all[1], WithinAbs(5.0, 1e-9));
    for (int32_t v = 0; v < prob.num_nodes(); ++v) {
        REQUIRE(all[v] >= first[v]);
    }
    const double depot_profit = prob.profit(0);
    const auto elim_first =
        cptp::preprocess::edge_elimination(prob, first, first, 0.0, depot_profit);
    const auto elim_all = cptp::preprocess::edge_elimination(prob, all, all, 0.0, depot_profit);
    constexpr int32_t kEdge12 = 2;
    REQUIRE_FALSE(elim_first[kEdge12]);
    REQUIRE(elim_all[kEdge12]);
}

TEST_CASE("Multicap: with one row the labeling is the same for both row options",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    const auto prob = make_problem(1);
    for (int32_t root = 0; root < prob.num_nodes(); ++root) {
        REQUIRE(cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kAll) ==
                cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kFirst));
    }
}

TEST_CASE("Multicap: a zero-demand negative cycle skips its row and terminates",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    // kRows[2] gives node 2 zero demand at profit 10: the walk 0-2-0 costs
    // 2 - 10 = -8 and adds no demand, so row 2's labeling never ends. With
    // no budget (0) it must still return, as row 2 skipped: rows 0 and 1.
    const auto prob = make_problem(3);
    const auto rows01 = make_problem(2);
    REQUIRE(cptp::preprocess::has_zero_demand_negative_cycle(prob, 2, prob.edge_costs(),
                                                             prob.profits()));
    REQUIRE_FALSE(cptp::preprocess::has_zero_demand_negative_cycle(prob, 0, prob.edge_costs(),
                                                                   prob.profits()));
    REQUIRE(cptp::preprocess::labeling_row(prob, 0, 2, prob.edge_costs(), prob.profits()).empty());
    const auto all = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kAll);
    REQUIRE_FALSE(all.empty());
    REQUIRE(all == cptp::preprocess::labeling_from(rows01, 0, 0, LabelRows::kAll));
}

TEST_CASE("Multicap: negative demands are rejected, a looping row 0 eliminates nothing",
          "[multicap][labeling]") {
    cptp::Problem bad;
    REQUIRE_THROWS_AS(bad.build(4, kEdges, kCosts, kProfits,
                                std::vector<std::vector<double>>{{0.0, 1.0, -1.0, 1.0}},
                                std::vector<double>{3.0}, 0, 0),
                      std::invalid_argument);
    // kRows[2] as row 0: the zero-demand walk 0-2-0 costs -8, so there are
    // no labeling bounds and the self-contained elimination fixes no edge.
    cptp::Problem loop;
    loop.build(4, kEdges, kCosts, kProfits, std::vector<std::vector<double>>{kRows[2]},
               std::vector<double>{kCaps[2]}, 0, 0);
    REQUIRE(cptp::preprocess::labeling_from(loop, 0).empty());
    const auto elim = cptp::preprocess::edge_elimination(loop, 0.0);
    REQUIRE(std::ranges::none_of(elim, [](bool b) { return b; }));
}

TEST_CASE("Multicap: Model and loaders reject a negative demand", "[multicap][io]") {
    // #95: with d = (0, 3, -2) and Q = 1 the heuristic took {1, 2} (load 1)
    // while the LP row 3 y_1 <= 1 forbade it. No path admits such a row.
    const auto negative = Catch::Matchers::Message("demands must be non-negative");
    cptp::Model model;
    model.set_graph(3, std::vector<cptp::Edge>{{0, 1}, {0, 2}, {1, 2}},
                    std::vector<double>{1.0, 1.0, 1.0});
    model.set_profits(std::vector<double>{0.0, 10.0, 10.0});
    model.set_depot(0);
    model.add_capacity_resource(std::vector<double>{0.0, 3.0, -2.0}, 1.0);
    REQUIRE_THROWS_MATCHES(model.solve(kQuiet), std::invalid_argument, negative);

    const auto tsplib = write_temp("cptp_multicap_negative.vrp",
                                   "NAME : negative\nTYPE : CVRP\nDIMENSION : 3\n"
                                   "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 1\n"
                                   "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                   "DEMAND_SECTION\n1 0\n2 3\n3 -2\n"
                                   "DEPOT_SECTION\n1\n-1\nEOF\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_tsplib(tsplib), std::invalid_argument, negative);
    std::filesystem::remove(tsplib);

    // The loader zeroes the depot's demand, but not a negative one.
    const auto depot = write_temp("cptp_multicap_negative_depot.vrp",
                                  "NAME : depot\nTYPE : CVRP\nDIMENSION : 3\n"
                                  "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 1\n"
                                  "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                  "DEMAND_SECTION\n1 -5\n2 1\n3 1\n"
                                  "DEPOT_SECTION\n1\n-1\nEOF\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_tsplib(depot), std::invalid_argument, negative);
    std::filesystem::remove(depot);

    // The numeric layout of the m = 2 test above, with d^1_2 = -1.
    const auto numeric = write_temp("cptp_multicap_negative.txt",
                                    "3 6 2\n"
                                    "0 0 0 0\n"
                                    "1 5 2 3\n"
                                    "2 4 1 -1\n"
                                    "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n"
                                    "4 7\n"
                                    "0 0\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(numeric), std::invalid_argument, negative);
    std::filesystem::remove(numeric);
}

TEST_CASE("Problem::build and the TSPLIB loader reject a source/target out of range", "[io]") {
    // #124: DEPOT_SECTION 5 with DIMENSION 3 stored source = target = 4.
    const auto range = Catch::Matchers::Message("source and target must be in [0, num_nodes)");
    const auto rows = std::vector<std::vector<double>>{kRows[0]};
    const auto caps = std::vector<double>{kCaps[0]};
    for (const auto& [s, t] :
         {std::pair{4, 4}, std::pair{0, 4}, std::pair{4, 0}, std::pair{-1, 0}, std::pair{0, -1}}) {
        INFO("source " << s << ", target " << t);
        cptp::Problem bad;
        REQUIRE_THROWS_MATCHES(bad.build(4, kEdges, kCosts, kProfits, rows, caps, s, t),
                               std::invalid_argument, range);
    }
    cptp::Problem path;
    path.build(4, kEdges, kCosts, kProfits, rows, caps, 0, 3);
    REQUIRE(path.target() == 3);

    const auto tsplib = write_temp("cptp_depot_out_of_range.vrp",
                                   "NAME : depot\nTYPE : CVRP\nDIMENSION : 3\n"
                                   "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 1\n"
                                   "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                   "DEMAND_SECTION\n1 0\n2 1\n3 1\n"
                                   "DEPOT_SECTION\n5\n-1\nEOF\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_tsplib(tsplib), std::invalid_argument, range);
    std::filesystem::remove(tsplib);
}

TEST_CASE("The TSPLIB loader rejects a DEPOT_SECTION id below 1", "[io]") {
    // #126: ids 0 and -3 were skipped, leaving the depot at node 0.
    const auto below = Catch::Matchers::Message("TSPLIB: DEPOT_SECTION id must be at least 1");
    for (const int id : {0, -3}) {
        INFO("depot id " << id);
        const auto path = write_temp("cptp_depot_below_one.vrp",
                                     "NAME : depot\nTYPE : CVRP\nDIMENSION : 3\n"
                                     "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 1\n"
                                     "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
                                     "DEMAND_SECTION\n1 0\n2 1\n3 1\n"
                                     "DEPOT_SECTION\n" +
                                         std::to_string(id) + "\n-1\nEOF\n");
        REQUIRE_THROWS_MATCHES(cptp::io::load_tsplib(path), std::runtime_error, below);
        std::filesystem::remove(path);
    }
}

TEST_CASE("The TSPLIB loader rejects a DEPOT_SECTION line that is not an int", "[io]") {
    // #129: an id that failed to parse was skipped, leaving the depot at node 0.
    const std::string head =
        "NAME : depot\nTYPE : CVRP\nDIMENSION : 3\n"
        "EDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 1\n"
        "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n"
        "DEMAND_SECTION\n1 0\n2 1\n3 1\n"
        "DEPOT_SECTION\n";
    for (const std::string id : {"abc", "99999999999", "1.5", "2x", "30 40"}) {
        INFO("depot line " << id);
        const auto path = write_temp("cptp_depot_not_int.vrp", head + id + "\n-1\nEOF\n");
        REQUIRE_THROWS_MATCHES(
            cptp::io::load_tsplib(path), std::runtime_error,
            Catch::Matchers::Message("TSPLIB: DEPOT_SECTION line \"" + id + "\" is not an int id"));
        std::filesystem::remove(path);
    }

    // The section ends at -1: an unknown section after it is not read as ids.
    const auto after = write_temp("cptp_depot_then_display.vrp",
                                  head + "2\n-1\nDISPLAY_DATA_SECTION\n1 0 0\n2 1 0\n3 0 1\nEOF\n");
    const auto prob = cptp::io::load_tsplib(after);
    std::filesystem::remove(after);
    REQUIRE(prob.source() == 1);
    REQUIRE(prob.target() == 1);
}

TEST_CASE("The TSPLIB loader sizes nothing by DIMENSION before its data", "[io]") {
    // #134: DIMENSION 200000 with no data allocated a 200000^2 matrix, and a
    // short FULL_MATRIX section was read out of bounds.
    const auto rejects = [](const std::string& name, const std::string& text,
                            const std::string& message) {
        INFO(name);
        const auto path = write_temp(name, text);
        REQUIRE_THROWS_MATCHES(cptp::io::load_tsplib(path), std::runtime_error,
                               Catch::Matchers::Message(message));
        std::filesystem::remove(path);
    };
    const std::string head = "NAME : t\nTYPE : CVRP\n";
    const std::string coords = "NODE_COORD_SECTION\n1 0 0\n2 1 0\n3 0 1\n";
    rejects("cptp_tsplib_huge_no_data.vrp", head + "DIMENSION : 200000\nEOF\n",
            "TSPLIB: no EDGE_WEIGHT_SECTION or NODE_COORD_SECTION");
    rejects("cptp_tsplib_huge_coords.vrp", head + "DIMENSION : 200000\n" + coords + "EOF\n",
            "TSPLIB: NODE_COORD_SECTION has 3 nodes, DIMENSION is 200000");
    rejects("cptp_tsplib_max_coords.vrp", head + "DIMENSION : 2147483647\n" + coords + "EOF\n",
            "TSPLIB: NODE_COORD_SECTION has 3 nodes, DIMENSION is 2147483647");
    rejects("cptp_tsplib_short_coords.vrp", head + "DIMENSION : 4\n" + coords + "EOF\n",
            "TSPLIB: NODE_COORD_SECTION has 3 nodes, DIMENSION is 4");

    const std::string full = head + "DIMENSION : 3\nEDGE_WEIGHT_TYPE : EXPLICIT\n";
    rejects("cptp_tsplib_short_full.vrp",
            full + "EDGE_WEIGHT_FORMAT : FULL_MATRIX\nEDGE_WEIGHT_SECTION\n0 1 2\n1 0 3\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 6 weights, FULL_MATRIX with DIMENSION 3 needs 9");
    rejects("cptp_tsplib_short_default_full.vrp", full + "EDGE_WEIGHT_SECTION\n0 1 2\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 3 weights, FULL_MATRIX with DIMENSION 3 needs 9");
    rejects("cptp_tsplib_huge_full.vrp",
            head + "DIMENSION : 2147483647\nEDGE_WEIGHT_SECTION\n0 1 2\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 3 weights, FULL_MATRIX with DIMENSION 2147483647 "
            "needs 4611686014132420609");
    rejects("cptp_tsplib_short_lower.vrp",
            full + "EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW\nEDGE_WEIGHT_SECTION\n0\n1 0\n2 3\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 5 weights, LOWER_DIAG_ROW with DIMENSION 3 needs 6");
    rejects("cptp_tsplib_upper_row.vrp",
            full + "EDGE_WEIGHT_FORMAT : UPPER_ROW\nEDGE_WEIGHT_SECTION\n1 2\n3\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_FORMAT UPPER_ROW is not supported");
    rejects("cptp_tsplib_short_node_weights.vrp",
            head + "DIMENSION : 3\n" + coords + "NODE_WEIGHT_SECTION\n0 -1\nEOF\n",
            "TSPLIB: NODE_WEIGHT_SECTION has 2 values, DIMENSION is 3");
    rejects("cptp_tsplib_short_profits.vrp",
            head + "DIMENSION : 3\n" + coords + "PROFIT_SECTION\n1 0\n2 5\nEOF\n",
            "TSPLIB: PROFIT_SECTION has 2 values, DIMENSION is 3");

    // An unknown section ends the one before it: its data does not make up
    // a short section, and a section with too many entries is an error too.
    const std::string display = "DISPLAY_DATA_SECTION\n1 0 0\n2 1 0\n3 0 1\n";
    rejects("cptp_tsplib_short_full_display.vrp",
            full + "EDGE_WEIGHT_SECTION\n0 1 2\n1 0 3\n" + display + "EOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 6 weights, FULL_MATRIX with DIMENSION 3 needs 9");
    rejects("cptp_tsplib_short_coords_display.vrp",
            head + "DIMENSION : 3\nNODE_COORD_SECTION\n1 0 0\n2 1 0\n" + display + "EOF\n",
            "TSPLIB: NODE_COORD_SECTION has 2 nodes, DIMENSION is 3");
    rejects("cptp_tsplib_long_full.vrp",
            full + "EDGE_WEIGHT_SECTION\n0 1 2\n1 0 3\n2 3 0\n4\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 10 weights, FULL_MATRIX with DIMENSION 3 needs 9");
    rejects("cptp_tsplib_long_coords.vrp", head + "DIMENSION : 2\n" + coords + "EOF\n",
            "TSPLIB: NODE_COORD_SECTION has 3 nodes, DIMENSION is 2");

    // A full section loads, and an unknown section after it is not read.
    const auto ok = write_temp("cptp_tsplib_full_ok.vrp",
                               full +
                                   "EDGE_WEIGHT_SECTION\n0 1 2\n1 0 3\n2 3 0\n"
                                   "DISPLAY_DATA_SECTION\n1 0 0\n2 1 0\n3 0 1\nEOF\n");
    const auto prob = cptp::io::load_tsplib(ok);
    std::filesystem::remove(ok);
    REQUIRE(prob.num_nodes() == 3);
    REQUIRE(prob.edge_costs() == std::vector<double>{1.0, 2.0, 3.0});
    const auto lower =
        write_temp("cptp_tsplib_lower_ok.vrp", full +
                                                   "EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW\n"
                                                   "EDGE_WEIGHT_SECTION\n0\n1 0\n2 3 0\nEOF\n");
    const auto lower_prob = cptp::io::load_tsplib(lower);
    std::filesystem::remove(lower);
    REQUIRE(lower_prob.edge_costs() == std::vector<double>{1.0, 2.0, 3.0});
    // LOWER_ROW has no diagonal: n(n-1)/2 weights (#142).
    const auto lower_row =
        write_temp("cptp_tsplib_lower_row_ok.vrp", full +
                                                       "EDGE_WEIGHT_FORMAT : LOWER_ROW\n"
                                                       "EDGE_WEIGHT_SECTION\n1\n2 3\nEOF\n");
    const auto lower_row_prob = cptp::io::load_tsplib(lower_row);
    std::filesystem::remove(lower_row);
    REQUIRE(lower_row_prob.edge_costs() == std::vector<double>{1.0, 2.0, 3.0});
    rejects("cptp_tsplib_short_lower_row.vrp",
            full + "EDGE_WEIGHT_FORMAT : LOWER_ROW\nEDGE_WEIGHT_SECTION\n1 2\nEOF\n",
            "TSPLIB: EDGE_WEIGHT_SECTION has 2 weights, LOWER_ROW with DIMENSION 3 needs 3");
}

TEST_CASE("The numeric loader rejects a source with no target", "[io]") {
    // #127: a lone trailing "2" loaded as an s-t path from node 2 to node 0.
    const std::string body =
        "3 6\n"
        "0 0 0\n1 5 2\n2 4 1\n"
        "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n"
        "4\n";
    const auto lone = write_temp("cptp_numeric_lone_source.txt", body + "2\n");
    REQUIRE_THROWS_MATCHES(
        cptp::io::load_numeric(lone), std::runtime_error,
        Catch::Matchers::Message("numeric: expected \"source target\", found a source only"));
    std::filesystem::remove(lone);

    const auto pair = write_temp("cptp_numeric_source_target.txt", body + "2 1\n");
    const auto path = cptp::io::load_numeric(pair);
    std::filesystem::remove(pair);
    REQUIRE(path.source() == 2);
    REQUIRE(path.target() == 1);

    const auto none = write_temp("cptp_numeric_no_source.txt", body);
    const auto tour = cptp::io::load_numeric(none);
    std::filesystem::remove(none);
    REQUIRE(tour.source() == 0);
    REQUIRE(tour.target() == 0);
}

TEST_CASE("The numeric loader rejects a capacity or source/target that does not parse", "[io]") {
    // #130: a failed read fell back to "no limit" or a tour from node 0.
    const std::string arcs =
        "3 6\n"
        "0 0 0\n1 5 2\n2 4 1\n"
        "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n";
    const auto rejects = [](const std::string& name, const std::string& text,
                            const std::string& message) {
        INFO(name);
        const auto path = write_temp(name, text);
        REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(path), std::runtime_error,
                               Catch::Matchers::Message(message));
        std::filesystem::remove(path);
    };
    const auto source = [](const std::string& t) {
        return "numeric: the source \"" + t + "\" is not an int";
    };
    const auto target = [](const std::string& t) {
        return "numeric: the target \"" + t + "\" is not an int";
    };
    const auto capacity = [](const std::string& t) {
        return "numeric: the capacity \"" + t + "\" is not a number";
    };
    rejects("cptp_numeric_text_source.txt", arcs + "4\nabc 1\n", source("abc"));
    rejects("cptp_numeric_overflow_source.txt", arcs + "4\n99999999999 1\n", source("99999999999"));
    // A failed >> consumes an overflowing or sign-only token: at end of file
    // these read as "no source line".
    rejects("cptp_numeric_overflow_source_end.txt", arcs + "4\n99999999999\n",
            source("99999999999"));
    rejects("cptp_numeric_sign_source_end.txt", arcs + "4\n-\n", source("-"));
    rejects("cptp_numeric_plus_minus_source.txt", arcs + "4\n+-1 1\n", source("+-1"));
    rejects("cptp_numeric_text_target.txt", arcs + "4\n2 x\n", target("x"));
    rejects("cptp_numeric_overflow_target_end.txt", arcs + "4\n2 99999999999\n",
            target("99999999999"));
    // With m = 1 a capacity token that is not a number hid the "s t" line.
    rejects("cptp_numeric_text_capacity.txt", arcs + "abc\n0 2\n", capacity("abc"));
    rejects("cptp_numeric_text_capacity_end.txt", arcs + "abc\n", capacity("abc"));
    rejects("cptp_numeric_overflow_capacity_end.txt", arcs + "1e999\n", capacity("1e999"));
    rejects("cptp_numeric_exponent_capacity_end.txt", arcs + "1e\n", capacity("1e"));
    rejects("cptp_numeric_sign_capacity_end.txt", arcs + "-\n", capacity("-"));

    // "+2" still reads as 2, as it did with operator>>.
    const auto plus = write_temp("cptp_numeric_plus_source.txt", arcs + "4\n+2 +1\n");
    const auto path = cptp::io::load_numeric(plus);
    std::filesystem::remove(plus);
    REQUIRE(path.source() == 2);
    REQUIRE(path.target() == 1);

    const auto frac = write_temp("cptp_numeric_fractional_capacity.txt", arcs + "2.5\n");
    const auto half = cptp::io::load_numeric(frac);
    std::filesystem::remove(frac);
    REQUIRE(half.capacities() == std::vector<double>{2.5});

    // No capacity line at all is still "no limit", and a tour.
    const auto none = write_temp("cptp_numeric_no_capacity.txt", arcs);
    const auto prob = cptp::io::load_numeric(none);
    std::filesystem::remove(none);
    REQUIRE(prob.capacities() == std::vector<double>{1e18});
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 0);
}

TEST_CASE("The numeric loader rejects a short arc list", "[io]") {
    // #128: num_arcs = 6 with 5 arc lines read an uninitialized head.
    const auto expected =
        Catch::Matchers::Message("numeric: expected 6 arc lines of \"tail head cost\"");
    const auto short_arcs = write_temp("cptp_numeric_short_arcs.txt",
                                       "3 6\n"
                                       "0 0 0\n1 5 2\n2 4 1\n"
                                       "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(short_arcs), std::runtime_error, expected);
    std::filesystem::remove(short_arcs);

    // A truncated last arc line, cut off before its cost.
    const auto cut = write_temp("cptp_numeric_cut_arc.txt",
                                "3 6\n"
                                "0 0 0\n1 5 2\n2 4 1\n"
                                "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(cut), std::runtime_error, expected);
    std::filesystem::remove(cut);

    // A short list followed by the capacity and "s t" lines: "10 0 2" was
    // read as the sixth arc and the file loaded as an unlimited tour.
    const auto trailing = write_temp("cptp_numeric_short_arcs_trailing.txt",
                                     "3 6\n"
                                     "0 0 0\n1 5 2\n2 4 1\n"
                                     "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n"
                                     "10\n0 2\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(trailing), std::runtime_error, expected);
    std::filesystem::remove(trailing);

    // An arc line with a fourth field.
    const auto extra = write_temp("cptp_numeric_extra_arc_field.txt",
                                  "3 6\n"
                                  "0 0 0\n1 5 2\n2 4 1\n"
                                  "0 1 1\n1 0 1 7\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n"
                                  "10\n");
    REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(extra), std::runtime_error, expected);
    std::filesystem::remove(extra);

    // Blank lines between arcs are still accepted.
    const auto blank = write_temp("cptp_numeric_blank_arc_lines.txt",
                                  "3 6\n"
                                  "0 0 0\n1 5 2\n2 4 1\n\n"
                                  "0 1 1\n1 0 1\n\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n"
                                  "10\n0 2\n");
    const auto prob = cptp::io::load_numeric(blank);
    std::filesystem::remove(blank);
    REQUIRE(prob.num_edges() == 3);
    REQUIRE(prob.capacities() == std::vector<double>{10.0});
    REQUIRE(prob.target() == 2);
}

TEST_CASE("The numeric loader rejects lines beyond the declared arcs", "[io]") {
    // #131: a sixth arc line with num_arcs = 5 loaded as capacity 2 and s-t
    // 1 -> 1, and the real "10" and "0 2" lines were ignored.
    const auto rejects = [](const std::string& name, const std::string& text,
                            const std::string& message) {
        INFO(name);
        const auto path = write_temp(name, text);
        REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(path), std::runtime_error,
                               Catch::Matchers::Message(message));
        std::filesystem::remove(path);
    };
    const std::string nodes = "0 0 0\n1 5 2\n2 4 1\n";
    const std::string five = "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n";
    rejects("cptp_numeric_extra_arc.txt", "3 5\n" + nodes + five + "2 1 1\n10\n0 2\n",
            "numeric: expected 1 capacities on one line, found 3");

    // m = 3: the extra arc reads as three capacities, and the real capacity
    // line as "s t".
    const std::string nodes3 = "0 0 0 0 0\n1 5 2 2 2\n2 4 1 1 1\n";
    rejects("cptp_numeric_extra_arc_m3.txt", "3 5 3\n" + nodes3 + five + "2 1 1\n10 10 10\n0 2\n",
            "numeric: expected \"source target\", found 3 fields");

    // Trailing fields on the "s t" line, and a line after it.
    const std::string body = "3 6\n" + nodes + five + "2 1 1\n10\n";
    rejects("cptp_numeric_st_extra_field.txt", body + "0 2 7\n",
            "numeric: expected \"source target\", found 3 fields");
    rejects("cptp_numeric_after_st.txt", body + "0 2\n\n5\n",
            R"(numeric: unexpected line after "source target": "5")");

    // Capacity and "s t" on one line, and m = 2 capacities split over lines.
    rejects("cptp_numeric_cap_st_one_line.txt", body.substr(0, body.size() - 3) + "10 0 2\n",
            "numeric: expected 1 capacities on one line, found 3");
    const std::string nodes2 = "0 0 0 0\n1 5 2 3\n2 4 1 6\n";
    rejects("cptp_numeric_split_caps.txt", "3 5 2\n" + nodes2 + five + "4\n7\n",
            "numeric: expected 2 capacities on one line, found 1");

    // Trailing blank lines are still accepted.
    const auto ok = write_temp("cptp_numeric_trailing_blank.txt", body + "0 2\n\n \t\n");
    const auto prob = cptp::io::load_numeric(ok);
    std::filesystem::remove(ok);
    REQUIRE(prob.capacities() == std::vector<double>{10.0});
    REQUIRE(prob.source() == 0);
    REQUIRE(prob.target() == 2);
}

TEST_CASE("The numeric loader rejects a bad m and node lines of the wrong width", "[io]") {
    // #132: "3 6 x" loaded with m = 1, "3 6 2.5" with m = 2, and an extra
    // node column shifted every later field.
    const auto rejects = [](const std::string& name, const std::string& text,
                            const std::string& message) {
        INFO(name);
        const auto path = write_temp(name, text);
        REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(path), std::runtime_error,
                               Catch::Matchers::Message(message));
        std::filesystem::remove(path);
    };
    const std::string nodes = "0 0 0\n1 5 2\n2 4 1\n";
    const std::string arcs = "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n";
    const auto m = [](const std::string& t) {
        return "numeric: the number of capacity rows \"" + t + "\" is not a positive int";
    };
    rejects("cptp_numeric_text_m.txt", "3 6 x\n" + nodes + arcs, m("x"));
    rejects("cptp_numeric_overflow_m.txt", "3 6 99999999999\n" + nodes + arcs, m("99999999999"));
    rejects("cptp_numeric_fractional_m.txt", "3 6 2.5\n" + nodes + arcs, m("2.5"));
    rejects("cptp_numeric_zero_m.txt", "3 6 0\n" + nodes + arcs, m("0"));
    rejects("cptp_numeric_negative_m.txt", "3 6 -1\n" + nodes + arcs, m("-1"));
    const std::string first = R"(numeric: first line must be "num_nodes num_arcs [m]")";
    rejects("cptp_numeric_first_extra.txt", "3 6 1 7\n" + nodes + arcs, first);
    rejects("cptp_numeric_first_short.txt", "3\n" + nodes + arcs, first);
    rejects("cptp_numeric_first_text.txt", "3x 6\n" + nodes + arcs, first);

    const std::string width = R"(numeric: expected 3 node lines of "id profit" and 1 demands)";
    rejects("cptp_numeric_node_extra.txt", "3 6\n0 0 0\n1 5 2 9\n2 4 1\n" + arcs, width);
    rejects("cptp_numeric_node_short.txt", "3 6\n0 0 0\n1 5\n2 4 1\n" + arcs, width);
    rejects("cptp_numeric_node_count.txt", "3 0\n0 0 0\n1 5 2\n", width);
    rejects("cptp_numeric_node_id.txt", "3 6\n0 0 0\nx 5 2\n2 4 1\n" + arcs,
            R"(numeric: the node id "x" is not an int)");
    rejects("cptp_numeric_node_demand.txt", "3 6\n0 0 0\n1 5 2x\n2 4 1\n" + arcs,
            R"(numeric: the node field "2x" is not a number)");

    // m = 2: two demands per node; blank lines between nodes still load.
    const auto ok = write_temp("cptp_numeric_m2_nodes.txt",
                               "3 6 2\n\n0 0 0 0\n1 5 2 3\n\n2 4 1 6\n" + arcs + "4 7\n");
    const auto two = cptp::io::load_numeric(ok);
    std::filesystem::remove(ok);
    REQUIRE(two.capacities() == std::vector<double>{4.0, 7.0});
    REQUIRE(two.demands(1) == std::vector<double>{0.0, 3.0, 6.0});
    REQUIRE(two.profits() == std::vector<double>{0.0, 5.0, 4.0});
}

TEST_CASE("The numeric loader rejects negative counts and sizes nothing by a count", "[io]") {
    // #133: "-1 6" threw length_error, "3 -1" loaded with no arcs, and a huge
    // m or num_nodes was allocated before the lines it counts were read.
    const auto rejects = [](const std::string& name, const std::string& text,
                            const std::string& message) {
        INFO(name);
        const auto path = write_temp(name, text);
        REQUIRE_THROWS_MATCHES(cptp::io::load_numeric(path), std::runtime_error,
                               Catch::Matchers::Message(message));
        std::filesystem::remove(path);
    };
    const std::string nodes = "0 0 0\n1 5 2\n2 4 1\n";
    const std::string arcs = "0 1 1\n1 0 1\n0 2 1\n2 0 1\n1 2 1\n2 1 1\n";
    const std::string negative = "numeric: num_nodes and num_arcs must not be negative";
    rejects("cptp_numeric_negative_nodes.txt", "-1 6\n" + nodes + arcs, negative);
    rejects("cptp_numeric_negative_arcs.txt", "3 -1\n" + nodes, negative);

    // The largest counts fail at the first line they count, not in an
    // allocation.
    rejects("cptp_numeric_huge_m.txt", "3 6 2147483647\n" + nodes + arcs,
            R"(numeric: expected 3 node lines of "id profit" and 2147483647 demands)");
    rejects("cptp_numeric_huge_m_no_nodes.txt", "0 0 2147483647\n",
            "numeric: expected 2147483647 capacities");
    rejects("cptp_numeric_huge_m_short_caps.txt", "0 0 2147483647\n4 7\n",
            "numeric: expected 2147483647 capacities on one line, found 2");
    rejects("cptp_numeric_huge_nodes.txt", "2147483647 6\n" + nodes + arcs,
            R"(numeric: expected 2147483647 node lines of "id profit" and 1 demands)");
    rejects("cptp_numeric_huge_arcs.txt", "3 2147483647\n" + nodes + arcs,
            R"(numeric: expected 2147483647 arc lines of "tail head cost")");
}

TEST_CASE("Problem::build rejects an edge endpoint out of range", "[io]") {
    // #125: edge (1, 5) with num_nodes 4 wrote past static_graph's degree array.
    const auto range = Catch::Matchers::Message("edge endpoints must be in [0, num_nodes)");
    const auto rows = std::vector<std::vector<double>>{kRows[0]};
    const auto caps = std::vector<double>{kCaps[0]};
    for (const auto& [u, v] :
         {std::pair{1, 5}, std::pair{5, 1}, std::pair{1, 4}, std::pair{-1, 2}, std::pair{2, -1}}) {
        INFO("edge (" << u << ", " << v << ")");
        auto edges = kEdges;
        edges.back() = {u, v};
        cptp::Problem bad;
        REQUIRE_THROWS_MATCHES(bad.build(4, edges, kCosts, kProfits, rows, caps, 0, 0),
                               std::invalid_argument, range);
    }
}

TEST_CASE("Multicap: zero-demand cycles that are not negative keep their row",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    // Row 1 gives customers 1..3 zero demand (the depot has 1), edges cost 3
    // and profits are 3: every zero-demand 2-cycle costs 6 - 6 = 0, so the
    // labeling terminates on dominance and row 1 is kept (f = 0 at 1..3,
    // where row 0 alone reaches 3 - 3 = 0 too).
    cptp::Problem prob;
    prob.build(4, kEdges, std::vector<double>(6, 3.0), std::vector<double>{0.0, 3.0, 3.0, 3.0},
               std::vector<std::vector<double>>{{0.0, 1.0, 1.0, 1.0}, {1.0, 0.0, 0.0, 0.0}},
               std::vector<double>{3.0, 1.0}, 0, 0);
    REQUIRE_FALSE(cptp::preprocess::has_zero_demand_negative_cycle(prob, 1, prob.edge_costs(),
                                                                   prob.profits()));
    const auto row1 = cptp::preprocess::labeling_row(prob, 0, 1, prob.edge_costs(), prob.profits());
    REQUIRE_FALSE(row1.empty());
    REQUIRE_THAT(row1[1], WithinAbs(0.0, 1e-9));
    REQUIRE_FALSE(cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kAll).empty());
}

TEST_CASE("Multicap: the surrogate row fixes an edge no single row does", "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    // Depot 0, customers 1, 2 (profit 1) and 3 (profit 15); edge (0,3) costs
    // 20, the others 1. Q = 10 on both rows, d0 = (0,4,7,5), d1 = (0,7,4,5):
    // 0-1-3 fits row 0 only (9 / 12), 0-2-3 row 1 only (12 / 9), so each row
    // alone reaches 3 at 2 - 16 = -14. The surrogate row d0/10 + d1/10 =
    // (0, 1.1, 1.1, 1.0) against 2 cuts both (2.1), leaving 0-3: 20 - 15 = 5.
    // No tour with 3 fits both rows, so the best is the empty tour (0) and
    // edge (1,3), lb = f[1] + 1 + f[3] = 0 + 1 + 5 = 6 > 0, is fixed.
    cptp::Problem prob;
    prob.build(4, kEdges, std::vector<double>{1.0, 1.0, 20.0, 1.0, 1.0, 1.0},
               std::vector<double>{0.0, 1.0, 1.0, 15.0},
               std::vector<std::vector<double>>{{0.0, 4.0, 7.0, 5.0}, {0.0, 7.0, 4.0, 5.0}},
               std::vector<double>{10.0, 10.0}, 0, 0);
    const auto all = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kAll);
    const auto sur = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kSurrogate);
    const auto first = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kFirst);
    REQUIRE_THAT(all[3], WithinAbs(-14.0, 1e-9));
    REQUIRE_THAT(sur[3], WithinAbs(5.0, 1e-9));
    for (int32_t v = 0; v < prob.num_nodes(); ++v) {
        REQUIRE(sur[v] >= first[v]);
    }
    constexpr int32_t kEdge13 = 4;
    const double depot_profit = prob.profit(0);
    REQUIRE_FALSE(cptp::preprocess::edge_elimination(prob, all, all, 0.0, depot_profit)[kEdge13]);
    REQUIRE(cptp::preprocess::edge_elimination(prob, sur, sur, 0.0, depot_profit)[kEdge13]);

    // Vector labels cut both routes too (each breaks one row): f[3] = 5.
    cptp::preprocess::LabelingStats stats;
    const auto vec = cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kVector, -1, &stats);
    REQUIRE_THAT(vec[3], WithinAbs(5.0, 1e-9));
    REQUIRE(stats.labels > 0);
    REQUIRE(cptp::preprocess::edge_elimination(prob, vec, vec, 0.0, depot_profit)[kEdge13]);
}

TEST_CASE("Multicap: vector labels are never weaker than per-row or surrogate",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    for (int32_t m = 1; m <= 3; ++m) {
        INFO("m = " << m);
        const auto prob = make_problem(m);
        for (int32_t root = 0; root < prob.num_nodes(); ++root) {
            const auto vec = cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kVector);
            const auto first = cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kFirst);
            if (m == 1) {
                REQUIRE(vec == first);
                continue;
            }
            const auto all = cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kAll);
            const auto sur = cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kSurrogate);
            for (int32_t v = 0; v < prob.num_nodes(); ++v) {
                REQUIRE(vec[v] >= all[v]);
                REQUIRE(vec[v] >= sur[v]);
            }
        }
    }
}

TEST_CASE("Multicap: the vector labeling gives no bound past its deadline (#115)",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    const auto prob = make_problem(3);
    const auto past = std::chrono::steady_clock::now();
    cptp::preprocess::LabelingStats stats;
    REQUIRE(cptp::preprocess::labeling_vector(prob, 0, prob.edge_costs(), prob.profits(), 0, &stats,
                                              past)
                .empty());
    REQUIRE(stats.timed_out);
    // (Row 0 stops at the same deadline, so labeling_from then gives no bound:
    // see the #143 test below.) A deadline that does not pass changes nothing.
    stats = {};
    const auto later = std::chrono::steady_clock::now() + std::chrono::hours(1);
    REQUIRE(cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kVector, -1, &stats, later) ==
            cptp::preprocess::labeling_from(prob, 0, 0, LabelRows::kVector));
    REQUIRE_FALSE(stats.timed_out);
}

TEST_CASE("Multicap: row 0's labeling gives no bound past its deadline (#143)",
          "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    const auto prob = make_problem(3);
    const auto past = std::chrono::steady_clock::now();
    bool timed_out = false;
    REQUIRE(cptp::preprocess::labeling_row(prob, 0, 0, prob.edge_costs(), prob.profits(), 0, past,
                                           &timed_out)
                .empty());
    REQUIRE(timed_out);
    // labeling_from then gives no bound either, under every row choice.
    for (const auto rows :
         {LabelRows::kFirst, LabelRows::kAll, LabelRows::kSurrogate, LabelRows::kVector}) {
        cptp::preprocess::LabelingStats stats;
        REQUIRE(cptp::preprocess::labeling_from(prob, 0, 0, rows, -1, &stats, past).empty());
        REQUIRE(stats.timed_out);
    }
    // A deadline that does not pass changes nothing.
    const auto later = std::chrono::steady_clock::now() + std::chrono::hours(1);
    for (const auto rows : {LabelRows::kFirst, LabelRows::kAll, LabelRows::kSurrogate}) {
        cptp::preprocess::LabelingStats stats;
        const auto bound = cptp::preprocess::labeling_from(prob, 0, 0, rows, -1, &stats, later);
        REQUIRE_FALSE(bound.empty());
        REQUIRE(bound == cptp::preprocess::labeling_from(prob, 0, 0, rows));
        REQUIRE_FALSE(stats.timed_out);
    }
}

TEST_CASE("Multicap: labeling_rows = vector keeps to time_limit (#115)",
          "[multicap][labeling][model]") {
    // The root vector labeling alone ran 972 s here under time_limit=120.
    constexpr double kLimit = 5.0;
    constexpr double kMargin = 5.0;
    auto prob = cptp::io::load("benchmarks/instances/multicap/spprclib/m10/B-n57-k7-20-m10.sppcc");
    cptp::Model model;
    model.set_problem(std::move(prob));
    const cptp::SolverOptions opts = {{"time_limit", std::to_string(kLimit)},
                                      {"output_flag", "false"},
                                      {"threads", "1"},
                                      {"labeling_rows", "vector"}};
    const auto start = std::chrono::steady_clock::now();
    (void)model.solve(opts);
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    REQUIRE(elapsed < kLimit + kMargin);
}

TEST_CASE("Multicap: the surrogate row needs two bounded rows", "[multicap][labeling]") {
    using cptp::preprocess::LabelRows;
    REQUIRE(cptp::preprocess::surrogate_row(make_problem(1)).demand.empty());
    const auto prob = make_problem(1);
    for (int32_t root = 0; root < prob.num_nodes(); ++root) {
        REQUIRE(cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kSurrogate) ==
                cptp::preprocess::labeling_from(prob, root, 0, LabelRows::kFirst));
    }
    // Rows 0 and 1 qualify (Q = 3, 3): s = d0/3 + d1/3 against 2 (+ slack).
    const auto row = cptp::preprocess::surrogate_row(make_problem(2));
    REQUIRE_THAT(row.demand[1], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(row.capacity, WithinAbs(2.0, 1e-6));
}

TEST_CASE("Multicap: labeling_rows = all keeps every row's optimum", "[multicap][labeling]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        INFO("m = " << c.rows);
        for (const char* rows : {"all", "surrogate", "vector"}) {
            INFO("labeling_rows = " << rows);
            cptp::Model model;
            fill_model(model, c.rows);
            model.set_depot(0);
            auto opts = kQuiet;
            opts.emplace_back("labeling_rows", rows);
            const auto result = model.solve(opts);
            REQUIRE(result.is_optimal());
            REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
            require_tour_fits(model.problem(), result.tour);
        }
    }
    REQUIRE_THROWS_AS(cptp::preprocess::parse_label_rows("every"), std::invalid_argument);
}

namespace {

// The exact minimum net cost of a (possibly non-elementary) walk from root to
// each node under one row, by DP over (accumulated demand, node): with integer
// demands >= 1 every step adds demand, so level d only reads lower levels.
std::vector<double> brute_force_labeling(const cptp::Problem& prob, int32_t root) {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    const auto n = static_cast<size_t>(prob.num_nodes());
    const auto q = static_cast<size_t>(prob.capacity());
    const auto& graph = prob.graph();
    std::vector<std::vector<double>> g(q + 1, std::vector<double>(n, kInf));
    g[static_cast<size_t>(prob.demand(root))][static_cast<size_t>(root)] = -prob.profit(root);
    std::vector<double> best(n, kInf);
    for (size_t d = 0; d <= q; ++d) {
        for (int32_t v = 0; v < prob.num_nodes(); ++v) {
            const auto dv = static_cast<size_t>(prob.demand(v));
            if (dv <= d) {
                for (auto e : graph.incident_edges(v)) {
                    const int32_t u = graph.other_endpoint(e, v);
                    g[d][static_cast<size_t>(v)] = std::min(
                        g[d][static_cast<size_t>(v)],
                        g[d - dv][static_cast<size_t>(u)] + prob.edge_cost(e) - prob.profit(v));
                }
            }
            best[static_cast<size_t>(v)] =
                std::min(best[static_cast<size_t>(v)], g[d][static_cast<size_t>(v)]);
        }
    }
    return best;
}

}  // namespace

TEST_CASE("Labeling: every live label is extended (#111)", "[labeling]") {
    // Random complete graphs with integer demands >= 1, against the exact DP.
    // Queueing (node, index) let erase_if shift a queued label's index onto a
    // newer one, so a live label went unextended and f came out too high.
    std::mt19937 rng(111);
    std::uniform_int_distribution<int> cost(1, 20);
    std::uniform_int_distribution<int> profit(0, 30);
    std::uniform_int_distribution<int> dem(1, 4);
    constexpr int32_t kNodes = 7;
    for (int trial = 0; trial < 300; ++trial) {
        std::vector<cptp::Edge> edges;
        std::vector<double> costs;
        for (int32_t u = 0; u < kNodes; ++u) {
            for (int32_t v = u + 1; v < kNodes; ++v) {
                edges.push_back({u, v});
                costs.push_back(cost(rng));
            }
        }
        std::vector<double> profits(kNodes);
        std::vector<double> demands(kNodes);
        for (int32_t v = 0; v < kNodes; ++v) {
            profits[static_cast<size_t>(v)] = profit(rng);
            demands[static_cast<size_t>(v)] = dem(rng);
        }
        cptp::Problem prob;
        prob.build(kNodes, edges, costs, profits, std::vector<std::vector<double>>{demands},
                   std::vector<double>{12.0}, 0, 0);
        const auto got =
            cptp::preprocess::labeling_row(prob, 0, 0, prob.edge_costs(), prob.profits());
        const auto want = brute_force_labeling(prob, 0);
        INFO("trial " << trial);
        REQUIRE(got.size() == want.size());
        for (size_t v = 0; v < want.size(); ++v) {
            INFO("node " << v);
            REQUIRE_THAT(got[v], WithinAbs(want[v], 1e-9));
        }
    }
}

namespace {

// The exact minimum net cost of a walk from root to each node at each load
// under one row, as brute_force_labeling, indexed g[d][v]. With two_cycle the
// walk never goes u -> v -> u: the DP then runs over (load, node, predecessor).
std::vector<std::vector<double>> brute_force_loads(const cptp::Problem& prob, int32_t root,
                                                   bool two_cycle) {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    const int32_t n = prob.num_nodes();
    const auto q = static_cast<size_t>(prob.capacity());
    const auto& graph = prob.graph();
    const auto none = static_cast<size_t>(n);  // the root's predecessor
    // h[d][v][p], p == none for the root label
    std::vector h(q + 1, std::vector(static_cast<size_t>(n), std::vector<double>(none + 1, kInf)));
    h[static_cast<size_t>(prob.demand(root))][static_cast<size_t>(root)][none] = -prob.profit(root);
    for (size_t d = 0; d <= q; ++d) {
        for (int32_t u = 0; u < n; ++u) {
            for (size_t p = 0; p <= none; ++p) {
                const double cost = h[d][static_cast<size_t>(u)][p];
                if (cost >= kInf) {
                    continue;
                }
                for (auto e : graph.incident_edges(u)) {
                    const int32_t v = graph.other_endpoint(e, u);
                    if (two_cycle && static_cast<size_t>(v) == p) {
                        continue;
                    }
                    const size_t dv = d + static_cast<size_t>(prob.demand(v));
                    if (dv > q) {
                        continue;
                    }
                    auto& next = h[dv][static_cast<size_t>(v)][static_cast<size_t>(u)];
                    next = std::min(next, cost + prob.edge_cost(e) - prob.profit(v));
                }
            }
        }
    }
    std::vector g(q + 1, std::vector<double>(static_cast<size_t>(n), kInf));
    for (size_t d = 0; d <= q; ++d) {
        for (size_t v = 0; v < static_cast<size_t>(n); ++v) {
            g[d][v] = *std::ranges::min_element(h[d][v]);
        }
    }
    return g;
}

// The exact load-consistent join of two such tables over edge e = (u, v):
// min over d1 + d2 <= cap of gf[d1][u] + c_e + gb[d2][v] + correction, both
// orientations.
double brute_force_join(const cptp::Problem& prob, int32_t e,
                        const std::vector<std::vector<double>>& gf,
                        const std::vector<std::vector<double>>& gb, double cap, double correction) {
    const auto& graph = prob.graph();
    const auto u = static_cast<size_t>(graph.edge_source(e));
    const auto v = static_cast<size_t>(graph.edge_target(e));
    double best = std::numeric_limits<double>::infinity();
    for (size_t d1 = 0; d1 < gf.size(); ++d1) {
        for (size_t d2 = 0; d2 < gb.size(); ++d2) {
            if (static_cast<double>(d1 + d2) > cap) {
                continue;
            }
            best = std::min({best, gf[d1][u] + gb[d2][v], gf[d1][v] + gb[d2][u]});
        }
    }
    return best + prob.edge_cost(e) + correction;
}

}  // namespace

TEST_CASE("Labeling: load join and 2-cycle elimination match an exact DP (#144)", "[labeling]") {
    // Random complete graphs with integer demands >= 1, as the #111 test, for
    // a tour (s = t = 0) and a path (s = 0, t = 6), with and without 2-cycle
    // elimination: f is the DP's per-node minimum and every edge's bound is
    // the exact join of the two halves' loads.
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::mt19937 rng(144);
    std::uniform_int_distribution<int> cost(1, 20);
    std::uniform_int_distribution<int> profit(0, 30);
    std::uniform_int_distribution<int> dem(1, 4);
    constexpr int32_t kNodes = 7;
    constexpr double kQ = 12.0;
    for (int trial = 0; trial < 200; ++trial) {
        std::vector<cptp::Edge> edges;
        std::vector<double> costs;
        for (int32_t u = 0; u < kNodes; ++u) {
            for (int32_t v = u + 1; v < kNodes; ++v) {
                edges.push_back({u, v});
                costs.push_back(cost(rng));
            }
        }
        std::vector<double> profits(kNodes);
        std::vector<double> demands(kNodes);
        for (int32_t v = 0; v < kNodes; ++v) {
            profits[static_cast<size_t>(v)] = profit(rng);
            demands[static_cast<size_t>(v)] = dem(rng);
        }
        for (const int32_t target : {0, kNodes - 1}) {
            cptp::Problem prob;
            prob.build(kNodes, edges, costs, profits, std::vector<std::vector<double>>{demands},
                       std::vector<double>{kQ}, 0, target);
            const bool tour = prob.is_tour();
            const double correction = tour ? prob.profit(0) : 0.0;
            const double cap = kQ + (tour ? prob.demand(0) : 0.0);
            for (const bool two_cycle : {false, true}) {
                INFO("trial " << trial << " target " << target << " two_cycle " << two_cycle);
                cptp::preprocess::LabelFronts ff;
                cptp::preprocess::LabelFronts fb;
                const auto f = cptp::preprocess::labeling_row(
                    prob, 0, 0, prob.edge_costs(), prob.profits(), 0,
                    cptp::preprocess::kNoLabelingDeadline, nullptr, two_cycle, &ff);
                const auto b =
                    tour ? f
                         : cptp::preprocess::labeling_row(
                               prob, target, 0, prob.edge_costs(), prob.profits(), 0,
                               cptp::preprocess::kNoLabelingDeadline, nullptr, two_cycle, &fb);
                const auto gf = brute_force_loads(prob, 0, two_cycle);
                const auto gb = tour ? gf : brute_force_loads(prob, target, two_cycle);
                for (size_t v = 0; v < static_cast<size_t>(kNodes); ++v) {
                    INFO("node " << v);
                    double want = kInf;
                    for (const auto& level : gf) {
                        want = std::min(want, level[v]);
                    }
                    REQUIRE_THAT(f[v], WithinAbs(want, 1e-9));
                }
                const auto lb = cptp::preprocess::edge_lower_bounds(prob, f, b, correction, &ff,
                                                                    tour ? &ff : &fb);
                for (auto e : prob.graph().edges()) {
                    INFO("edge " << e);
                    const double want = brute_force_join(prob, e, gf, gb, cap, correction);
                    if (want >= kInf) {
                        REQUIRE(lb[static_cast<size_t>(e)] >= kInf);
                    } else {
                        REQUIRE_THAT(lb[static_cast<size_t>(e)], WithinAbs(want, 1e-9));
                    }
                }
            }
        }
    }
}

TEST_CASE("Multicap: labeling_join = load and labeling_2cycle keep the optimum (#144)",
          "[multicap][labeling]") {
    struct Case {
        int32_t rows;
        double objective;
    };
    for (const auto& c : {Case{1, -29.0}, Case{2, -20.0}, Case{3, -18.0}}) {
        INFO("m = " << c.rows);
        for (const char* join : {"min", "load"}) {
            for (const char* two_cycle : {"false", "true"}) {
                INFO("labeling_join = " << join << ", labeling_2cycle = " << two_cycle);
                cptp::Model model;
                fill_model(model, c.rows);
                model.set_depot(0);
                auto opts = kQuiet;
                opts.emplace_back("labeling_join", join);
                opts.emplace_back("labeling_2cycle", two_cycle);
                const auto result = model.solve(opts);
                REQUIRE(result.is_optimal());
                REQUIRE_THAT(result.objective, WithinAbs(c.objective, 1e-6));
                require_tour_fits(model.problem(), result.tour);
            }
        }
    }
    cptp::Model model;
    fill_model(model, 1);
    model.set_depot(0);
    auto opts = kQuiet;
    opts.emplace_back("labeling_join", "max");
    REQUIRE_THROWS_AS(model.solve(opts), std::invalid_argument);
}

TEST_CASE("Labeling join: s-t path solves agree across labeling_join and 2cycle (#144)",
          "[multicap][labeling]") {
    // Random s-t path instances (s = 0, t = 7) with a binding capacity: the
    // load join and 2-cycle elimination must not cut off the optimum that the
    // plain join (labeling_join = min) finds.
    std::mt19937 rng(1440);
    std::uniform_int_distribution<int> cost(1, 20);
    std::uniform_int_distribution<int> profit(0, 30);
    std::uniform_int_distribution<int> dem(1, 4);
    constexpr int32_t kNodes = 8;
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<cptp::Edge> edges;
        std::vector<double> costs;
        for (int32_t u = 0; u < kNodes; ++u) {
            for (int32_t v = u + 1; v < kNodes; ++v) {
                edges.push_back({u, v});
                costs.push_back(cost(rng));
            }
        }
        std::vector<double> profits(kNodes);
        std::vector<double> demands(kNodes);
        for (int32_t v = 0; v < kNodes; ++v) {
            profits[static_cast<size_t>(v)] = profit(rng);
            demands[static_cast<size_t>(v)] = dem(rng);
        }
        cptp::Problem prob;
        prob.build(kNodes, edges, costs, profits, std::vector<std::vector<double>>{demands},
                   std::vector<double>{10.0}, 0, kNodes - 1);
        auto solve = [&](const char* join, const char* two_cycle) {
            cptp::Model model;
            model.set_problem(prob);
            auto opts = kQuiet;
            opts.emplace_back("labeling_join", join);
            opts.emplace_back("labeling_2cycle", two_cycle);
            const auto result = model.solve(opts);
            REQUIRE(result.is_optimal());
            return result.objective;
        };
        INFO("trial " << trial);
        const double want = solve("min", "false");
        REQUIRE_THAT(solve("load", "false"), WithinAbs(want, 1e-6));
        REQUIRE_THAT(solve("min", "true"), WithinAbs(want, 1e-6));
        REQUIRE_THAT(solve("load", "true"), WithinAbs(want, 1e-6));
    }
}
