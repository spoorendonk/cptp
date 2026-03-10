#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/io.h"
#include "model/model.h"

using Catch::Matchers::WithinAbs;

static const cptp::SolverOptions quiet = {
    {"output_flag", "false"},
};

static cptp::SolveResult solve_instance(const char* path,
                                        int time_limit = 120) {
  auto prob = cptp::io::load(path);
  cptp::Model model;
  model.set_problem(std::move(prob));
  auto opts = quiet;
  opts.push_back({"time_limit", std::to_string(time_limit)});
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
  auto r =
      solve_instance("benchmarks/instances/spprclib/B-n45-k6-54.sppcc", 60);
  REQUIRE(r.is_optimal());
  REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
}

TEST_CASE("B-n50-k8-40 optimal", "[optimal][slow]") {
  auto r =
      solve_instance("benchmarks/instances/spprclib/B-n50-k8-40.sppcc", 120);
  REQUIRE(r.is_optimal());
  REQUIRE_THAT(r.objective, WithinAbs(-12832.0, 1.0));
}

TEST_CASE("B-n52-k7-15 optimal", "[optimal][slow]") {
  auto r =
      solve_instance("benchmarks/instances/spprclib/B-n52-k7-15.sppcc", 60);
  REQUIRE(r.is_optimal());
  REQUIRE_THAT(r.objective, WithinAbs(-74998.0, 1.0));
}

TEST_CASE("A-n63-k10-44 optimal", "[optimal][slow]") {
  auto r =
      solve_instance("benchmarks/instances/spprclib/A-n63-k10-44.sppcc", 60);
  REQUIRE(r.is_optimal());
  REQUIRE_THAT(r.objective, WithinAbs(-32561.0, 1.0));
}

TEST_CASE("A-n69-k9-42 optimal", "[optimal][slow]") {
  auto r =
      solve_instance("benchmarks/instances/spprclib/A-n69-k9-42.sppcc", 60);
  REQUIRE(r.is_optimal());
  REQUIRE_THAT(r.objective, WithinAbs(-43290.0, 1.0));
}
