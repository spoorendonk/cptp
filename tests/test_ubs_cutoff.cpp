#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/io.h"
#include "model/model.h"

using Catch::Matchers::WithinAbs;

static const cptp::SolverOptions quiet = {
    {"output_flag", "false"},
};

/// Solve with known UB cutoff and all heuristics disabled (prove-only mode).
static cptp::SolveResult prove_instance(const char* path, double cutoff,
                                          int time_limit = 120) {
    auto prob = cptp::io::load(path);
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", std::to_string(time_limit)});
    opts.push_back({"cutoff", std::to_string(cutoff)});
    opts.push_back({"heu_ws", "false"});
    return model.solve(opts);
}

/// Solve normally (baseline).
static cptp::SolveResult solve_instance(const char* path, int time_limit = 120) {
    auto prob = cptp::io::load(path);
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", std::to_string(time_limit)});
    return model.solve(opts);
}

// ── Tiny instances: prove-only should still find optimality ──

TEST_CASE("prove-only: tiny4 tour", "[prove]") {
    auto r = prove_instance("tests/data/tiny4.txt", -11.0, 30);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("prove-only: tiny4 path", "[prove]") {
    auto r = prove_instance("tests/data/tiny4_path.txt", -13.0, 30);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

// ── SPPRCLIB instances: prove-only with known optimal UBS ──

TEST_CASE("prove-only: B-n45-k6-54", "[prove][slow]") {
    auto r = prove_instance("benchmarks/instances/spprclib/B-n45-k6-54.sppcc", -74278.0, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

TEST_CASE("prove-only: B-n50-k8-40", "[prove][slow]") {
    auto r = prove_instance("benchmarks/instances/spprclib/B-n50-k8-40.sppcc", -12832.0, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-12832.0, 1.0));
}

TEST_CASE("prove-only: B-n52-k7-15", "[prove][slow]") {
    auto r = prove_instance("benchmarks/instances/spprclib/B-n52-k7-15.sppcc", -74998.0, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74998.0, 1.0));
}

TEST_CASE("prove-only: A-n63-k10-44", "[prove][slow]") {
    auto r = prove_instance("benchmarks/instances/spprclib/A-n63-k10-44.sppcc", -32561.0, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-32561.0, 1.0));
}

TEST_CASE("prove-only: A-n69-k9-42", "[prove][slow]") {
    auto r = prove_instance("benchmarks/instances/spprclib/A-n69-k9-42.sppcc", -43290.0, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-43290.0, 1.0));
}

// ── Comparison: prove-only vs normal should reach same optimum ──

TEST_CASE("prove-only matches normal solve: tiny4", "[prove]") {
    auto normal = solve_instance("tests/data/tiny4.txt", 30);
    REQUIRE(normal.is_optimal());

    auto prove = prove_instance("tests/data/tiny4.txt", normal.objective, 30);
    REQUIRE(prove.is_optimal());
    REQUIRE_THAT(prove.objective, WithinAbs(normal.objective, 1e-6));
}

// ── RC fixing tests ──

/// Solve with RC fixing enabled in prove-only mode.
static cptp::SolveResult prove_rc_fixing(const char* path, double cutoff,
                                           const char* strategy = "root_only",
                                           bool fix_to_one = false,
                                           int time_limit = 30) {
    auto prob = cptp::io::load(path);
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", std::to_string(time_limit)});
    opts.push_back({"cutoff", std::to_string(cutoff)});
    opts.push_back({"heu_ws", "false"});
    opts.push_back({"rc_fixing", strategy});
    if (fix_to_one)
        opts.push_back({"rc_fixing_to_one", "true"});
    return model.solve(opts);
}

TEST_CASE("rc_fixing root_only: tiny4 tour optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4.txt", -11.0);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing root_only: tiny4 path optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4_path.txt", -13.0);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

TEST_CASE("rc_fixing on_ub_improvement: tiny4 tour optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "on_ub_improvement");
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing with fix_to_one: tiny4 tour optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "root_only", true);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing matches baseline objective: tiny4", "[prove][rc_fixing]") {
    auto baseline = prove_instance("tests/data/tiny4.txt", -11.0, 30);
    REQUIRE(baseline.is_optimal());

    auto rc = prove_rc_fixing("tests/data/tiny4.txt", -11.0);
    REQUIRE(rc.is_optimal());
    REQUIRE_THAT(rc.objective, WithinAbs(baseline.objective, 1e-6));
}

TEST_CASE("rc_fixing periodic: tiny4 tour optimal", "[prove][rc_fixing]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", "30"});
    opts.push_back({"cutoff", "-11.0"});
    opts.push_back({"heu_ws", "false"});
    opts.push_back({"rc_fixing", "periodic"});
    opts.push_back({"rc_fixing_interval", "10"});
    auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing with fix_to_one: tiny4 path optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4_path.txt", -13.0, "root_only", true);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

TEST_CASE("rc_fixing on_ub_improvement: tiny4 path optimal", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4_path.txt", -13.0, "on_ub_improvement");
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

TEST_CASE("rc_fixing off matches on_ub_improvement on tiny4", "[prove][rc_fixing]") {
    auto r_off = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "off");
    auto r_default = prove_instance("tests/data/tiny4.txt", -11.0, 30);  // default = on_ub_improvement
    REQUIRE(r_off.is_optimal());
    REQUIRE(r_default.is_optimal());
    REQUIRE_THAT(r_off.objective, WithinAbs(r_default.objective, 1e-6));
}

TEST_CASE("rc_fixing all strategies agree on optimum: tiny4", "[prove][rc_fixing]") {
    auto r_root = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "root_only");
    auto r_ub   = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "on_ub_improvement");
    auto r_per  = prove_rc_fixing("tests/data/tiny4.txt", -11.0, "periodic");
    REQUIRE(r_root.is_optimal());
    REQUIRE(r_ub.is_optimal());
    REQUIRE(r_per.is_optimal());
    REQUIRE_THAT(r_root.objective, WithinAbs(r_ub.objective, 1e-6));
    REQUIRE_THAT(r_root.objective, WithinAbs(r_per.objective, 1e-6));
}

TEST_CASE("rc_fixing default active in normal solve: tiny4", "[rc_fixing]") {
    // Default is on_ub_improvement — should work without prove-only mode
    auto r = solve_instance("tests/data/tiny4.txt", 30);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing default active in normal solve: tiny4 path", "[rc_fixing]") {
    auto prob = cptp::io::load("tests/data/tiny4_path.txt");
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto r = model.solve(quiet);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

TEST_CASE("rc_fixing periodic: tiny4 path optimal", "[prove][rc_fixing]") {
    auto prob = cptp::io::load("tests/data/tiny4_path.txt");
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", "30"});
    opts.push_back({"cutoff", "-13.0"});
    opts.push_back({"heu_ws", "false"});
    opts.push_back({"rc_fixing", "periodic"});
    opts.push_back({"rc_fixing_interval", "10"});
    auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

TEST_CASE("rc_fixing periodic + fix_to_one: tiny4 tour", "[prove][rc_fixing]") {
    auto prob = cptp::io::load("tests/data/tiny4.txt");
    cptp::Model model;
    model.set_problem(std::move(prob));
    auto opts = quiet;
    opts.push_back({"time_limit", "30"});
    opts.push_back({"cutoff", "-11.0"});
    opts.push_back({"heu_ws", "false"});
    opts.push_back({"rc_fixing", "periodic"});
    opts.push_back({"rc_fixing_interval", "10"});
    opts.push_back({"rc_fixing_to_one", "true"});
    auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-11.0, 1.0));
}

TEST_CASE("rc_fixing on_ub_improvement + fix_to_one: tiny4 path", "[prove][rc_fixing]") {
    auto r = prove_rc_fixing("tests/data/tiny4_path.txt", -13.0,
                              "on_ub_improvement", true);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-13.0, 1.0));
}

// ── RC fixing on SPPRCLIB (slow) ──

TEST_CASE("rc_fixing root_only: B-n45-k6-54", "[prove][rc_fixing][slow]") {
    auto r = prove_rc_fixing(
        "benchmarks/instances/spprclib/B-n45-k6-54.sppcc", -74278.0,
        "root_only", false, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

TEST_CASE("rc_fixing on_ub_improvement: B-n45-k6-54", "[prove][rc_fixing][slow]") {
    auto r = prove_rc_fixing(
        "benchmarks/instances/spprclib/B-n45-k6-54.sppcc", -74278.0,
        "on_ub_improvement", false, 60);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}
