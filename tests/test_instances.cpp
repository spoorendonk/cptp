// End-to-end branch-and-cut benchmarks on SPPRCLIB instances.
//
// Loads real instances with built-in profits, solves them with the full B&C
// solver, and reports:
//   - Objective value, bound, optimality gap
//   - Total solve time, B&C nodes explored
//   - Cuts added per separator type, separation rounds
//
// Usage:
//   ./test_instances                         # all benchmarks
//   ./test_instances "[small]"               # small instances only
//   ./test_instances "[spprclib]"            # SPPRCLIB instances only
//   ./test_instances --benchmark-samples 3   # fewer samples for timing

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "core/io.h"
#include "core/problem.h"
#include "core/solution.h"
#include "model/model.h"

using namespace cptp;

namespace {

std::filesystem::path find_instances_dir() {
    for (const auto& candidate : {
        std::filesystem::path("benchmarks/instances"),
        std::filesystem::path("../benchmarks/instances"),
        std::filesystem::path("../../benchmarks/instances"),
    }) {
        if (std::filesystem::exists(candidate)) return candidate;
    }
    return "benchmarks/instances";
}

SolveResult solve_instance(Problem& prob, double time_limit) {
    Model model;
    model.set_problem(std::move(prob));

    SolverOptions opts = {
        {"time_limit", std::to_string(time_limit)},
        {"output_flag", "false"},
    };

    return model.solve(opts);
}

void print_result(const std::string& instance_name, const Problem& prob,
                  const SolveResult& result) {
    std::cout << "\n=== " << instance_name << " ===\n";
    std::cout << "  Nodes: " << prob.num_nodes()
              << "  Edges: " << prob.num_edges()
              << "  Capacity: " << prob.capacity() << "\n";

    const char* status_str = "???";
    switch (result.status) {
        case SolveResult::Status::Optimal:   status_str = "OPTIMAL"; break;
        case SolveResult::Status::Feasible:  status_str = "FEASIBLE"; break;
        case SolveResult::Status::TimeLimit: status_str = "TIME_LIMIT"; break;
        case SolveResult::Status::Infeasible: status_str = "INFEASIBLE"; break;
        case SolveResult::Status::Unbounded: status_str = "UNBOUNDED"; break;
        case SolveResult::Status::Error:     status_str = "ERROR"; break;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Status: " << status_str << "\n";
    std::cout << "  Objective: " << result.objective
              << "  Bound: " << result.bound
              << "  Gap: " << (result.gap * 100.0) << "%\n";
    std::cout << "  Time: " << result.time_seconds << "s"
              << "  B&C nodes: " << result.nodes << "\n";

    if (result.total_cuts > 0) {
        std::cout << "  Separation: " << result.total_cuts << " cuts in "
                  << result.separation_rounds << " rounds\n";
        std::cout << std::fixed << std::setprecision(3);
        for (const auto& [name, stats] : result.separator_stats) {
            std::cout << "    " << std::setw(10) << std::left << name
                      << std::right
                      << std::setw(6) << stats.cuts_added << " cuts"
                      << std::setw(6) << stats.rounds_called << " rounds"
                      << std::setw(8) << stats.time_seconds << "s\n";
        }
    }

    if (!result.tour.empty()) {
        std::cout << "  Tour (" << result.tour.size() - 1 << " stops): ";
        for (size_t i = 0; i < std::min(result.tour.size(), size_t(10)); ++i) {
            if (i > 0) std::cout << " -> ";
            std::cout << result.tour[i];
        }
        if (result.tour.size() > 10) std::cout << " -> ...";
        std::cout << "\n";
    }
}

}  // namespace

// ========== Small SPPRCLIB instances (< 55 nodes) ==========

TEST_CASE("E2E: B-n45-k6-54 (45 nodes)", "[small][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "B-n45-k6-54.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 60.0);
    print_result("B-n45-k6-54", prob, result);

    REQUIRE(result.has_solution());
    CHECK(result.total_cuts > 0);

    BENCHMARK("B-n45-k6-54") {
        auto p = io::load(path);
        return solve_instance(p, 60.0);
    };
}

TEST_CASE("E2E: P-n50-k7-92 (50 nodes)", "[small][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "P-n50-k7-92.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 60.0);
    print_result("P-n50-k7-92", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("P-n50-k7-92") {
        auto p = io::load(path);
        return solve_instance(p, 60.0);
    };
}

TEST_CASE("E2E: B-n52-k7-15 (52 nodes)", "[small][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "B-n52-k7-15.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 60.0);
    print_result("B-n52-k7-15", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("B-n52-k7-15") {
        auto p = io::load(path);
        return solve_instance(p, 60.0);
    };
}

TEST_CASE("E2E: A-n54-k7-149 (54 nodes)", "[small][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "A-n54-k7-149.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 60.0);
    print_result("A-n54-k7-149", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("A-n54-k7-149") {
        auto p = io::load(path);
        return solve_instance(p, 60.0);
    };
}

// ========== Medium SPPRCLIB instances (55-80 nodes) ==========

TEST_CASE("E2E: P-n55-k7-116 (55 nodes)", "[medium][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "P-n55-k7-116.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 120.0);
    print_result("P-n55-k7-116", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("P-n55-k7-116") {
        auto p = io::load(path);
        return solve_instance(p, 120.0);
    };
}

TEST_CASE("E2E: A-n63-k9-157 (63 nodes)", "[medium][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "A-n63-k9-157.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 120.0);
    print_result("A-n63-k9-157", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("A-n63-k9-157") {
        auto p = io::load(path);
        return solve_instance(p, 120.0);
    };
}

TEST_CASE("E2E: B-n68-k9-65 (68 nodes)", "[medium][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "B-n68-k9-65.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 120.0);
    print_result("B-n68-k9-65", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("B-n68-k9-65") {
        auto p = io::load(path);
        return solve_instance(p, 120.0);
    };
}

TEST_CASE("E2E: E-n76-k7-44 (76 nodes)", "[medium][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "E-n76-k7-44.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 120.0);
    print_result("E-n76-k7-44", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("E-n76-k7-44") {
        auto p = io::load(path);
        return solve_instance(p, 120.0);
    };
}

// ========== Large SPPRCLIB instances (100+ nodes) ==========

TEST_CASE("E2E: E-n101-k8-291 (101 nodes)", "[large][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "E-n101-k8-291.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 300.0);
    print_result("E-n101-k8-291", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("E-n101-k8-291") {
        auto p = io::load(path);
        return solve_instance(p, 300.0);
    };
}

TEST_CASE("E2E: M-n101-k10-97 (101 nodes)", "[large][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "M-n101-k10-97.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 300.0);
    print_result("M-n101-k10-97", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("M-n101-k10-97") {
        auto p = io::load(path);
        return solve_instance(p, 300.0);
    };
}

TEST_CASE("E2E: M-n121-k7-260 (121 nodes)", "[large][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "M-n121-k7-260.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 300.0);
    print_result("M-n121-k7-260", prob, result);

    REQUIRE(result.has_solution());

    BENCHMARK("M-n121-k7-260") {
        auto p = io::load(path);
        return solve_instance(p, 300.0);
    };
}

// ========== Extra-large (200+ nodes, long timeout) ==========

TEST_CASE("E2E: M-n200-k16-143 (200 nodes)", "[xlarge][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "M-n200-k16-143.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 600.0);
    print_result("M-n200-k16-143", prob, result);

    REQUIRE(result.has_solution());
}

TEST_CASE("E2E: G-n262-k25-316 (262 nodes)", "[xlarge][spprclib][e2e]") {
    auto path = find_instances_dir() / "spprclib" / "G-n262-k25-316.sppcc";
    REQUIRE(std::filesystem::exists(path));

    auto prob = io::load(path);
    auto result = solve_instance(prob, 600.0);
    print_result("G-n262-k25-316", prob, result);

    REQUIRE(result.has_solution());
}

