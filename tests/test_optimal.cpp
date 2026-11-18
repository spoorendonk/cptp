#include "core/io.h"
#include "model/model.h"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <regex>
#include <span>
#include <sstream>
#include <string>
#include <thread>

using Catch::Matchers::WithinAbs;

// The [optimal] tests solve with a threads option of at most 4, and CMake
// declares 4 PROCESSORS for them, so ctest -j runs fewer of them side by side
// and a time-limited solve stays clear of its limit (#160). The default,
// (hw+1)/2, grows with the host; below 8 logical CPUs it stays in force.
static const std::string kOptimalThreads =
    std::to_string(std::min(4U, (std::max(1U, std::thread::hardware_concurrency()) + 1U) / 2U));

static const cptp::SolverOptions kQuiet = {
    {"output_flag", "false"},
    {"threads", kOptimalThreads},
};

static cptp::SolveResult solve_instance(const char* path, int time_limit = 120) {
    auto prob = cptp::io::load(path);
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = kQuiet;
    opts.emplace_back("time_limit", std::to_string(time_limit));
    return model.solve(opts);
}

// --- Known optimal values ---
// Format: {instance_path, expected_objective}

TEST_CASE("tiny4 optimal", "[optimal]") {
    auto r = solve_instance("tests/data/tiny4.txt", 30);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("tiny4_path optimal", "[optimal][path]") {
    auto r = solve_instance("tests/data/tiny4_path.txt", 30);
    REQUIRE(r.is_optimal());
    // Path 0 -> 1 -> 3: cost = 10 + 7 = 17, profit = 0 + 20 + 10 = 30, obj = -13
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
    // Path should start at source=0 and end at target=3
    REQUIRE(r.tour.front() == 0);
    REQUIRE(r.tour.back() == 3);
}

TEST_CASE("B-n45-k6-54 optimal", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/B-n45-k6-54.sppcc", 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

TEST_CASE("B-n50-k8-40 optimal", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/B-n50-k8-40.sppcc", 120);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-12832.0, 1.0));
}

TEST_CASE("B-n52-k7-15 optimal", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/B-n52-k7-15.sppcc", 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74998.0, 1.0));
}

TEST_CASE("A-n63-k10-44 optimal", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/A-n63-k10-44.sppcc", 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-32561.0, 1.0));
}

TEST_CASE("A-n69-k9-42 optimal", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/A-n69-k9-42.sppcc", 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-43290.0, 1.0));
}

// Regression guard for the HiGHS 1.15 warm-start reroute: feeding the
// warm-start via setSolution() together with the LP-guided user-solution
// callback corrupted branch-and-bound on this instance and proved a wrong
// optimum (-30950 instead of -32628). Runs the default heu_ws + heu_lpg path,
// at the threads value kOptimalThreads gives (the default on hosts below 8
// logical CPUs).
TEST_CASE("M-n101-k10-97 optimal (warm-start reroute)", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/M-n101-k10-97.sppcc", 120);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-32628.0, 1.0));
}

namespace {

struct OptimumCase {
    const char* path;
    double objective;
};

// Eight instances whose optima the parallel tree search must reach.
constexpr std::array kParallelCases = {
    OptimumCase{"tests/data/tiny4.txt", -11.0},
    OptimumCase{"tests/data/tiny4_path.txt", -13.0},
    OptimumCase{"benchmarks/instances/spprclib/B-n45-k6-54.sppcc", -74278.0},
    OptimumCase{"benchmarks/instances/spprclib/B-n50-k8-40.sppcc", -12832.0},
    OptimumCase{"benchmarks/instances/spprclib/B-n52-k7-15.sppcc", -74998.0},
    OptimumCase{"benchmarks/instances/spprclib/A-n63-k10-44.sppcc", -32561.0},
    OptimumCase{"benchmarks/instances/spprclib/A-n69-k9-42.sppcc", -43290.0},
    OptimumCase{"benchmarks/instances/spprclib/M-n101-k10-97.sppcc", -32628.0},
};

}  // namespace

// Solves `cases` with HiGHS's parallel tree search on (plus `extra`) and
// checks each optimum, and that the main MIP ran more than one tree worker
// (a sub-MIP reports 1), so the parallel search really ran.
static void check_optima_under_parallel_tree_search(
    const cptp::SolverOptions& extra, std::span<const OptimumCase> cases = kParallelCases) {
    for (const auto& [path, objective] : cases) {
        INFO(path);
        auto prob = cptp::io::load(path);
        cptp::Model model;
        model.set_problem(std::move(prob));
        cptp::SolverOptions opts = {
            {"output_flag", "true"}, {"time_limit", "180"}, {"parallel", "on"}};
        opts.insert(opts.end(), extra.begin(), extra.end());
        // The Model's logger writes to std::cout, HiGHS's log included.
        std::ostringstream log;
        cptp::SolveResult r;
        {
            struct CoutRedirect {
                std::streambuf* old;
                ~CoutRedirect() { std::cout.rdbuf(old); }
            } redirect{std::cout.rdbuf(log.rdbuf())};
            r = model.solve(opts);
        }
        REQUIRE(r.is_optimal());
        REQUIRE_THAT(r.objective, WithinAbs(objective, 1.0));
        const std::string text = log.str();
        const std::regex workers(R"(using up to (\d+) workers across)");
        bool several_workers = false;
        for (auto it = std::sregex_iterator(text.begin(), text.end(), workers);
             it != std::sregex_iterator(); ++it) {
            several_workers = several_workers || std::stoi((*it)[1].str()) > 1;
        }
        REQUIRE(several_workers);
    }
}

// HiGHS's parallel tree search with its workers run one after another on the
// calling thread (mip_search_simulate_concurrency): every worker path of the
// patched callbacks runs, deterministically. A worker's incumbent must pass the
// lazy-constraint check before it tightens the worker's bounds, and an integral
// node whose LP solution has subtours must stay open (#152).
TEST_CASE("optimal under simulated parallel tree search", "[optimal][slow][parallel]") {
    check_optima_under_parallel_tree_search(
        {{"threads", "4"}, {"mip_search_simulate_concurrency", "true"}});
}

// Real parallel tree search: the workers call the bridge's callbacks at once
// (#153). Timing-dependent, so the gate solves one instance; the hidden tests
// below cover the rest.
TEST_CASE("optimal under parallel tree search", "[optimal][slow][parallel]") {
    check_optima_under_parallel_tree_search({{"threads", "4"}},
                                            std::span(kParallelCases).subspan(2, 1));
}

// Every instance at threads 4, hidden from the gate: the TSan stress run,
// ./cptp_tests "[parallel-stress]".
TEST_CASE("optimal under parallel tree search, all instances", "[.parallel-stress]") {
    check_optima_under_parallel_tree_search({{"threads", "4"}});
}

// The #153 acceptance matrix, hidden from the default run (and so from the
// gate): ./build/cptp_tests "[parallel-matrix]".
TEST_CASE("optimal under parallel tree search, threads x seeds", "[.parallel-matrix]") {
    const int threads = GENERATE(2, 4, 8);
    const int seed = GENERATE(0, 1, 2);
    INFO("threads=" << threads << " random_seed=" << seed);
    check_optima_under_parallel_tree_search(
        {{"threads", std::to_string(threads)}, {"random_seed", std::to_string(seed)}});
}

// Hyperplane branching through auxiliary integer columns (#137): every mode
// adds its columns and reaches the branch_hyper=off optimum.
TEST_CASE("B-n45-k6-54 optimal with hyperplane branching", "[optimal][slow][hyperplane]") {
    for (const char* mode : {"pairs", "clusters", "demand", "cardinality", "all"}) {
        INFO("branch_hyper=" << mode);
        auto prob = cptp::io::load("benchmarks/instances/spprclib/B-n45-k6-54.sppcc");
        cptp::Model model;
        model.set_problem(std::move(prob));
        cptp::SolverOptions opts = {{"output_flag", "false"},
                                    {"time_limit", "60"},
                                    {"threads", "1"},
                                    {"branch_hyper", mode}};
        auto r = model.solve(opts);
        REQUIRE(r.hyperplane_columns > 0);
        REQUIRE(r.is_optimal());
        REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
    }
}

// The benchmark pipeline reads per-phase timings and fixing counts from
// SolveResult; verify they are wired and internally consistent. Runs with
// fixing + propagation enabled so the propagator/rc-fixing paths are active.
TEST_CASE("SolveResult per-phase stats are populated and consistent", "[optimal][slow]") {
    auto prob = cptp::io::load("benchmarks/instances/spprclib/B-n45-k6-54.sppcc");
    cptp::Model model;
    model.set_problem(std::move(prob));
    cptp::SolverOptions opts = {{"output_flag", "false"},
                                {"time_limit", "60"},
                                {"threads", kOptimalThreads},
                                {"rc_fixing", "adaptive"},
                                {"bounds_propagation", "true"}};
    auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    // Wall-clock phase timers are non-negative.
    REQUIRE(r.separation_time_seconds >= 0.0);
    REQUIRE(r.propagator_time_seconds >= 0.0);
    REQUIRE(r.heuristic_time_seconds >= 0.0);
    REQUIRE(r.rc_time_seconds >= 0.0);
    // Separation runs on every CPTP solve, so its timer must have advanced.
    REQUIRE(r.separation_time_seconds > 0.0);
    // rc-fixing runs inside the propagator, so its time is a subset of it.
    REQUIRE(r.rc_time_seconds <= r.propagator_time_seconds + 1e-9);
    // Fixing counters are counts (non-negative).
    REQUIRE(r.sweep_fixings >= 0);
    REQUIRE(r.chain_fixings >= 0);
    REQUIRE(r.rc_fix0_count >= 0);
    REQUIRE(r.rc_fix1_count >= 0);
}
