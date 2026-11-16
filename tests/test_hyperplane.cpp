// Hyperplane branching through auxiliary integer columns (#137): one column
// z_h = a_h . y per hyperplane, which HiGHS branches on as on any column.
#include "core/problem.h"
#include "Highs.h"
#include "model/highs_bridge.h"
#include "model/hyperplane.h"
#include "model/model.h"
#include "util/logger.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using Catch::Matchers::WithinAbs;

namespace {

// Tour from depot 0 over customers 1..3. Pairs: 1's nearest customer is 2
// (cost 6), and 2 and 3 pair with an already paired node, so one pair (1, 2).
// Clusters: seeded by 2 (demand 4), grown to {2, 3, 1}, so one cluster.
const std::vector<cptp::Edge> kEdges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
const std::vector<double> kCosts = {10.0, 8.0, 12.0, 6.0, 7.0, 5.0};
const std::vector<double> kProfits = {0.0, 20.0, 15.0, 10.0};

cptp::SolveResult solve(const std::vector<double>& demands, double capacity,
                        const std::string& mode) {
    cptp::Model model;
    model.set_graph(4, kEdges, kCosts);
    model.set_depot(0);
    model.set_profits(kProfits);
    model.add_capacity_resource(demands, capacity);
    const cptp::SolverOptions opts = {
        {"output_flag", "false"}, {"time_limit", "30"}, {"threads", "1"}, {"branch_hyper", mode}};
    return model.solve(opts);
}

struct ModeColumns {
    const char* mode;
    int32_t columns;
};

}  // namespace

TEST_CASE("Hyperplane columns: each mode adds its columns and keeps the off optimum",
          "[hyperplane][branching]") {
    const std::vector<double> demands = {0.0, 3.0, 4.0, 2.0};
    const auto off = solve(demands, 7.0, "off");
    REQUIRE(off.is_optimal());
    REQUIRE(off.hyperplane_columns == 0);
    for (const auto [mode, columns] :
         {ModeColumns{"pairs", 1}, ModeColumns{"clusters", 1}, ModeColumns{"demand", 1},
          ModeColumns{"cardinality", 1}, ModeColumns{"all", 4}}) {
        INFO("branch_hyper=" << mode);
        const auto r = solve(demands, 7.0, mode);
        REQUIRE(r.hyperplane_columns == columns);
        REQUIRE(r.is_optimal());
        REQUIRE_THAT(r.objective, WithinAbs(off.objective, 1e-6));
    }
}

TEST_CASE("Hyperplane columns: non-integer demands give no demand-weighted column",
          "[hyperplane][branching]") {
    // z = 2.5 y_1 + ... is not integer at an integer y, so an integer z column
    // would cut off feasible tours. The row is skipped; pairs and cardinality,
    // whose coefficients are 1, stay.
    const std::vector<double> demands = {0.0, 2.5, 3.5, 1.5};
    const auto off = solve(demands, 6.0, "off");
    REQUIRE(off.is_optimal());
    for (const auto [mode, columns] :
         {ModeColumns{"demand", 0}, ModeColumns{"clusters", 0}, ModeColumns{"all", 2}}) {
        INFO("branch_hyper=" << mode);
        const auto r = solve(demands, 6.0, mode);
        REQUIRE(r.hyperplane_columns == columns);
        REQUIRE(r.is_optimal());
        REQUIRE_THAT(r.objective, WithinAbs(off.objective, 1e-6));
    }

    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, demands, 6.0, 0, 0);
    const cptp::detail::HyperplaneCandidates gen(prob, 6, "all");
    REQUIRE(gen.rows().empty());
    REQUIRE(gen.skipped_rows() == std::vector<int32_t>{0});
}

TEST_CASE("Hyperplane columns: a completed warm start satisfies every row",
          "[hyperplane][branching]") {
    cptp::Problem prob;
    prob.build(4, kEdges, kCosts, kProfits, std::vector<double>{0.0, 3.0, 4.0, 2.0}, 10.0, 0, 0);
    Highs highs;
    highs.setOptionValue("output_flag", false);
    cptp::Logger logger;
    cptp::HiGHSBridge bridge(prob, highs, logger);
    bridge.set_hyperplanes(
        cptp::detail::HyperplaneCandidates(prob, bridge.y_offset(), "all").columns());
    bridge.build_formulation();
    REQUIRE(bridge.num_hyperplane_columns() == 4);

    // Tour 0-1-2-0 over edges 0 (0,1), 1 (0,2) and 3 (1,2); y_0 = y_1 = y_2 = 1.
    const auto m = static_cast<int32_t>(kEdges.size());
    std::vector<double> sol(static_cast<size_t>(m + 4), 0.0);
    for (const int32_t e : {0, 1, 3}) {
        sol[static_cast<size_t>(e)] = 1.0;
    }
    for (const int32_t i : {0, 1, 2}) {
        sol[static_cast<size_t>(m) + static_cast<size_t>(i)] = 1.0;
    }
    bridge.complete_solution(sol);

    const HighsLp& lp = highs.getLp();
    REQUIRE(static_cast<HighsInt>(sol.size()) == lp.num_col_);
    // pair_1_2 = 2, cluster_0 = 3 + 4, demand = 3 + 4, cardinality = 2.
    const std::vector<double> z(sol.begin() + m + 4, sol.end());
    REQUIRE(z == std::vector<double>{2.0, 7.0, 7.0, 2.0});

    HighsSparseMatrix a = lp.a_matrix_;
    a.ensureColwise();  // addRow may leave the LP's matrix row-wise
    std::vector<double> activity(static_cast<size_t>(lp.num_row_), 0.0);
    for (HighsInt c = 0; c < lp.num_col_; ++c) {
        REQUIRE(sol[c] >= lp.col_lower_[c] - 1e-9);
        REQUIRE(sol[c] <= lp.col_upper_[c] + 1e-9);
        for (HighsInt k = a.start_[c]; k < a.start_[c + 1]; ++k) {
            activity[a.index_[k]] += a.value_[k] * sol[c];
        }
    }
    for (HighsInt r = 0; r < lp.num_row_; ++r) {
        INFO("row " << r);
        REQUIRE(activity[r] >= lp.row_lower_[r] - 1e-9);
        REQUIRE(activity[r] <= lp.row_upper_[r] + 1e-9);
    }
}
