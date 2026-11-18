#include "core/io.h"
#include "Highs.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsUserSeparator.h"
#include "model/model.h"
#include "model/options.h"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

// HiGHS's root cut-loop limits through the patch's hooks (#150).

namespace {

constexpr int kHighs = HighsUserHooks::kRootLimitHighs;
constexpr int kUnlimited = HighsUserHooks::kRootLimitUnlimited;

// The separator adds at most this many cuts, one per call.
constexpr int64_t kMaxCuts = 60;

struct RootRun {
    int64_t calls = 0;          // separator calls in the main MIP's root cut loop
    HighsInt seen_rounds = -2;  // rootSepaRounds(mipsolver, 7) inside the solve
    HighsInt seen_stall = -2;   // rootSepaStall(mipsolver, 3) inside the solve
    HighsStatus status = HighsStatus::kError;
};

// max x + y over x, y in [0, 2] with an integer z that nothing constrains.
// Call k of the separator adds x + (1 + k/1000) y <= 2 - k/1000, which the
// LP solution of the call before (x = 2 - (k-1)/1000, y = 0) violates by
// 1/1000. Every LP solution is integral; the lazy check rejects each one
// until all kMaxCuts cuts are in, so the root LP never closes the solve.
// The root cut loop runs the separator once per round and stops at its
// round cap, or after call kMaxCuts, whose LP solution then passes the
// check and closes the gap. HiGHS records
// rootlpsol right after the loop, so an empty one marks a call in it.
RootRun run_root_loop(int max_rounds, int stall_limit) {
    HighsLp lp;
    lp.num_col_ = 3;
    lp.num_row_ = 0;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = {1.0, 1.0, 0.0};
    lp.col_lower_ = {0.0, 0.0, 0.0};
    lp.col_upper_ = {2.0, 2.0, 1.0};
    lp.integrality_ = {HighsVarType::kContinuous, HighsVarType::kContinuous,
                       HighsVarType::kInteger};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 3;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0, 0, 0};

    RootRun run;
    int64_t cuts = 0;
    HighsUserSeparator::setCallback([&run, &cuts](const HighsLpRelaxation&, HighsCutPool& cutpool,
                                                  const HighsMipSolver& mipsolver) {
        if (mipsolver.submip) {
            return;
        }
        if (mipsolver.mipdata_->rootlpsol.empty()) {
            ++run.calls;
            run.seen_rounds = HighsUserSeparator::rootSepaRounds(mipsolver, 7);
            run.seen_stall = HighsUserSeparator::rootSepaStall(mipsolver, 3);
        }
        if (cuts == kMaxCuts) {
            return;
        }
        const auto k = static_cast<double>(++cuts);
        std::array<HighsInt, 2> idx = {0, 1};
        std::array<double, 2> val = {1.0, 1.0 + (k / 1000.0)};
        cutpool.addCut(mipsolver, idx.data(), val.data(), 2, 2.0 - (k / 1000.0));
    });
    HighsUserSeparator::setFeasibilityCheck(
        [&cuts](const std::vector<double>&) { return cuts == kMaxCuts; });
    HighsUserSeparator::setRootSeparationLimits(max_rounds, stall_limit);

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("time_limit", 60.0);
    // An earlier solve on this thread may have left a scheduler of another
    // size, and HiGHS refuses threads = 1 against it.
    Highs::resetThreadScheduler(true);
    if (highs.passModel(lp) == HighsStatus::kOk) {
        run.status = highs.run();
    }
    // Before any REQUIRE: the hooks capture this frame's run, and a later
    // test on this thread must not run them.
    HighsUserSeparator::clearCallback();
    return run;
}

}  // namespace

TEST_CASE("Root separation: limits are pending until clearCallback resets them",
          "[highs_patch][root_separation]") {
    HighsUserSeparator::clearCallback();
    REQUIRE(HighsUserHooks::capture() == nullptr);

    // Limits alone make a solve capture the hooks.
    HighsUserSeparator::setRootSeparationLimits(5, kUnlimited);
    const auto hooks = HighsUserHooks::capture();
    // Before any REQUIRE, so a failure leaves no limits pending.
    HighsUserSeparator::clearCallback();
    REQUIRE(hooks != nullptr);
    CHECK(hooks->root_sepa_rounds == 5);
    CHECK(hooks->root_sepa_stall == kUnlimited);

    CHECK(HighsUserHooks::pending().root_sepa_rounds == kHighs);
    CHECK(HighsUserHooks::pending().root_sepa_stall == kHighs);
    REQUIRE(HighsUserHooks::capture() == nullptr);
}

TEST_CASE("Root separation: the round cap reaches HiGHS's root cut loop",
          "[highs_patch][root_separation]") {
    SECTION("a cap of N runs N rounds") {
        for (const int cap : {1, 4, 9}) {
            INFO("cap " << cap);
            const RootRun run = run_root_loop(cap, kUnlimited);
            REQUIRE(run.status == HighsStatus::kOk);
            CHECK(run.seen_rounds == cap);
            CHECK(run.seen_stall == kHighsIInf);
            CHECK(run.calls == cap);
        }
    }
    SECTION("unlimited runs until the last cut is in") {
        const RootRun run = run_root_loop(kUnlimited, kUnlimited);
        REQUIRE(run.status == HighsStatus::kOk);
        CHECK(run.seen_rounds == kHighsIInf);
        CHECK(run.calls == kMaxCuts);
    }
    SECTION("highs keeps HiGHS's cap, below the cuts on offer") {
        const RootRun run = run_root_loop(kHighs, kHighs);
        REQUIRE(run.status == HighsStatus::kOk);
        CHECK(run.seen_rounds == 7);
        CHECK(run.seen_stall == 3);
        CHECK(run.calls < kMaxCuts);
    }
    SECTION("a stall limit reaches the solve") {
        const RootRun run = run_root_loop(kHighs, 11);
        REQUIRE(run.status == HighsStatus::kOk);
        CHECK(run.seen_rounds == 7);
        CHECK(run.seen_stall == 11);
    }
}

TEST_CASE("Root separation: the stall limit ends HiGHS's root cut loop",
          "[highs_patch][root_separation]") {
    // Each cut moves the bound by 1/1000, which HiGHS counts as a stall once
    // the first rounds are past: from there each extra stalled round allowed
    // is one more separator call, until the last cut is in.
    const RootRun highs3 = run_root_loop(kUnlimited, kHighs);
    REQUIRE(highs3.status == HighsStatus::kOk);
    CHECK(highs3.seen_stall == 3);
    CHECK(highs3.calls < kMaxCuts);
    CHECK(run_root_loop(kUnlimited, 3).calls == highs3.calls);
    for (const int stall : {5, 11}) {
        INFO("stall " << stall);
        const RootRun run = run_root_loop(kUnlimited, stall);
        REQUIRE(run.status == HighsStatus::kOk);
        CHECK(run.calls == highs3.calls + (stall - 3));
    }
    CHECK(run_root_loop(kUnlimited, kUnlimited).calls == kMaxCuts);
}

// The options reach HiGHS through Model::solve: on a root-only solve
// (mip_max_nodes 1, optimum -74998 as the cutoff), a round cap of 1 leaves a
// weaker root bound than HiGHS's cap, and no cap closes the root (Mac mini:
// -219062, -88009 and -74998, each in under a second).
TEST_CASE("Root separation: Model::solve hands the options to HiGHS (#150)",
          "[model][root_separation]") {
    auto root_bound = [](const char* rounds) {
        cptp::Model model;
        model.set_problem(cptp::io::load("benchmarks/instances/spprclib/B-n52-k7-15.sppcc"));
        const cptp::SolverOptions opts = {
            {"output_flag", "false"},
            {"time_limit", "60"},
            {"threads", "1"},
            {"mip_max_nodes", "1"},
            {"heu_ws", "false"},
            {"heu_lpg", "false"},
            {"cutoff", "-74997.9"},
            {"root_sepa_rounds", rounds},
            {"root_sepa_stall", "unlimited"},
        };
        const auto r = model.solve(opts);
        REQUIRE(r.nodes <= 1);
        return r.bound;
    };
    const double one = root_bound("1");
    const double highs = root_bound("highs");
    const double unlimited = root_bound("unlimited");
    INFO("bounds: 1 " << one << ", highs " << highs << ", unlimited " << unlimited);
    CHECK(one < highs - 1.0);
    CHECK(highs < unlimited - 1.0);
    CHECK(unlimited > -74998.5);
}

TEST_CASE("Root separation: root_sepa_rounds and root_sepa_stall parse (#150)",
          "[model][options][root_separation]") {
    for (const char* key : {"root_sepa_rounds", "root_sepa_stall"}) {
        INFO(key);
        CHECK(cptp::parse_root_limit_option(key, "highs") == kHighs);
        CHECK(cptp::parse_root_limit_option(key, "unlimited") == kUnlimited);
        CHECK(cptp::parse_root_limit_option(key, "25") == 25);
        for (const char* v : {"0", "-1", "", "Highs", "inf", "3x"}) {
            INFO(v);
            CHECK_THROWS_AS(cptp::parse_root_limit_option(key, v), std::invalid_argument);
        }
    }
    CHECK(cptp::root_limit_name(kHighs) == "highs");
    CHECK(cptp::root_limit_name(kUnlimited) == "unlimited");
    CHECK(cptp::root_limit_name(25) == "25");

    cptp::Model model;
    const std::vector<cptp::Edge> edges = {{0, 1}, {0, 2}, {1, 2}};
    model.set_graph(3, edges, std::vector<double>{5.0, 3.0, 4.0});
    model.set_depot(0);
    model.set_profits(std::vector<double>{0.0, 10.0, 8.0});
    model.add_capacity_resource(std::vector<double>{0.0, 2.0, 3.0}, 5.0);
    const cptp::SolverOptions quiet = {{"output_flag", "false"}, {"time_limit", "30"}};
    for (const char* v : {"highs", "unlimited", "3"}) {
        INFO(v);
        auto opts = quiet;
        opts.emplace_back("root_sepa_rounds", v);
        opts.emplace_back("root_sepa_stall", v);
        CHECK(model.solve(opts).is_optimal());
    }
    auto bad = quiet;
    bad.emplace_back("root_sepa_rounds", "0");
    CHECK_THROWS_AS(model.solve(bad), std::invalid_argument);
}
