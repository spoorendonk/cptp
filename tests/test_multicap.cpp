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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
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
// x12 = x23 = 1 (each customer has degree 2).
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

LpPoint make_point(const cptp::Problem& prob) {
    LpPoint p;
    const auto& g = prob.graph();
    p.x.assign(static_cast<size_t>(prob.num_edges()), 0.0);
    cptp::digraph_builder builder(prob.num_nodes());
    for (auto e : g.edges()) {
        const int32_t u = g.edge_source(e);
        const int32_t v = g.edge_target(e);
        const double xe = point_value(std::min(u, v), std::max(u, v));
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

std::vector<cptp::sep::Cut> cuts_of(const cptp::Problem& prob, cptp::sep::Separator& separator) {
    const auto point = make_point(prob);
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

// The candidates of `mode` at y = 0.5 everywhere (y columns from offset 6).
constexpr int32_t kYOffset = 6;

std::vector<cptp::detail::Hyperplane> hyperplanes_at_half(const cptp::Problem& prob,
                                                          const std::string& mode) {
    const cptp::detail::HyperplaneCandidates gen(prob, kYOffset, mode);
    const std::vector<double> sol(static_cast<size_t>(kYOffset + prob.num_nodes()), 0.5);
    return gen(sol);
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
    const auto hs = hyperplanes_at_half(prob, "all");
    const auto names = names_of(hs);
    // One cluster {1, 2, 3}: row 0's load 1.5 and row 1's 2.5 are fractional.
    for (const char* name : {"cluster_0", "cluster_0_row1", "demand", "demand_row1"}) {
        INFO(name);
        REQUIRE(std::ranges::count(names, std::string(name)) == 1);
    }
    const auto row1 =
        std::ranges::find(hs, std::string("demand_row1"), &cptp::detail::Hyperplane::name);
    REQUIRE(row1 != hs.end());
    REQUIRE(row1->indices == std::vector<int32_t>{7, 8, 9});
    REQUIRE(row1->values == std::vector<double>{2.0, 2.0, 1.0});
}

TEST_CASE("Multicap: one row keeps the one-row hyperplane names", "[multicap][branching]") {
    const auto prob = make_problem(1);
    const auto demand = hyperplanes_at_half(prob, "demand");
    REQUIRE(names_of(demand) == std::vector<std::string>{"demand"});
    REQUIRE(demand[0].indices == std::vector<int32_t>{7, 8, 9});
    REQUIRE(demand[0].values == std::vector<double>{1.0, 1.0, 1.0});
    REQUIRE(names_of(hyperplanes_at_half(prob, "clusters")) ==
            std::vector<std::string>{"cluster_0"});
    for (const auto& name : names_of(hyperplanes_at_half(prob, "all"))) {
        REQUIRE_FALSE(name.contains("_row"));
    }
}

TEST_CASE("Multicap: rows with Q_k = 0 or no limit give no hyperplane", "[multicap][branching]") {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, {kRows[0], kRows[1], kRows[2], kRows[1]},
               {3.0, 0.0, cptp::sep::kInfiniteCapacity, 3.0}, 0, 0);
    const cptp::detail::HyperplaneCandidates gen(prob, kYOffset, "demand");
    REQUIRE(gen.rows() == std::vector<int32_t>{0, 3});
    REQUIRE(names_of(hyperplanes_at_half(prob, "demand")) ==
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
        for (const char* rows : {"all", "surrogate"}) {
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
