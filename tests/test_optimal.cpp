#include "core/io.h"
#include "model/model.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

using Catch::Matchers::WithinAbs;

static const cptp::SolverOptions kQuiet = {
    {"output_flag", "false"},
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
// optimum (-30950 instead of -32628). Runs the default heu_ws + heu_lpg path.
TEST_CASE("M-n101-k10-97 optimal (warm-start reroute)", "[optimal][slow]") {
    auto r = solve_instance("benchmarks/instances/spprclib/M-n101-k10-97.sppcc", 120);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-32628.0, 1.0));
}

// The bridge's callback state is not yet safe under concurrent tree workers
// (#153), so Model::solve keeps HiGHS's tree search serial even when asked for
// HiGHS's parallel MIP (the default from v1.16).
TEST_CASE("B-n45-k6-54 optimal with parallel requested", "[optimal][slow]") {
    auto prob = cptp::io::load("benchmarks/instances/spprclib/B-n45-k6-54.sppcc");
    cptp::Model model;
    model.set_problem(std::move(prob));
    cptp::SolverOptions opts = {
        {"output_flag", "false"}, {"time_limit", "60"}, {"parallel", "on"}, {"threads", "4"}};
    auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

// HiGHS's parallel tree search with its workers run one after another on the
// calling thread (mip_search_simulate_concurrency): every worker path of the
// patched callbacks runs, deterministically. A worker's incumbent must pass the
// lazy-constraint check before it tightens the worker's bounds, and an integral
// node whose LP solution has subtours must stay open (#152).
TEST_CASE("optimal under simulated parallel tree search", "[optimal][slow][parallel]") {
    struct Case {
        const char* path;
        double objective;
    };
    for (const auto& [path, objective] : {
             Case{"tests/data/tiny4.txt", -11.0},
             Case{"tests/data/tiny4_path.txt", -13.0},
             Case{"benchmarks/instances/spprclib/B-n45-k6-54.sppcc", -74278.0},
             Case{"benchmarks/instances/spprclib/B-n50-k8-40.sppcc", -12832.0},
             Case{"benchmarks/instances/spprclib/B-n52-k7-15.sppcc", -74998.0},
             Case{"benchmarks/instances/spprclib/A-n63-k10-44.sppcc", -32561.0},
             Case{"benchmarks/instances/spprclib/A-n69-k9-42.sppcc", -43290.0},
             Case{"benchmarks/instances/spprclib/M-n101-k10-97.sppcc", -32628.0},
         }) {
        INFO(path);
        auto prob = cptp::io::load(path);
        cptp::Model model;
        model.set_problem(std::move(prob));
        const cptp::SolverOptions opts = {{"output_flag", "true"},
                                          {"time_limit", "180"},
                                          {"parallel", "on"},
                                          {"threads", "4"},
                                          {"mip_search_simulate_concurrency", "true"}};
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
        // The main MIP runs more than one tree worker (a sub-MIP reports 1), so
        // the parallel option was not forced off.
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
