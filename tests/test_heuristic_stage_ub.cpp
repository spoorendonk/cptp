#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "core/io.h"
#include "heuristic/primal_heuristic.h"

namespace {

std::vector<std::filesystem::path> collect_heuristic_gate_instances() {
  std::vector<std::filesystem::path> paths;
  auto collect = [&](const std::filesystem::path& dir, const char* ext) {
    if (!std::filesystem::exists(dir)) return;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      if (!entry.is_regular_file()) continue;
      if (entry.path().extension() == ext) {
        paths.push_back(entry.path());
      }
    }
  };

  collect("benchmarks/instances/spprclib", ".sppcc");
  collect("benchmarks/instances/roberti", ".vrp");
  std::sort(paths.begin(), paths.end());
  return paths;
}

std::filesystem::path find_snapshot_csv() {
  for (const auto& candidate : {
           std::filesystem::path("tests/data/heuristic_stage_ub_ils.csv"),
           std::filesystem::path("../tests/data/heuristic_stage_ub_ils.csv"),
           std::filesystem::path("../../tests/data/heuristic_stage_ub_ils.csv"),
       }) {
    if (std::filesystem::exists(candidate)) return candidate;
  }
  return "tests/data/heuristic_stage_ub_ils.csv";
}

std::map<std::string, double> load_expected_snapshot() {
  const auto csv_path = find_snapshot_csv();
  std::ifstream in(csv_path);
  if (!in) return {};

  std::map<std::string, double> snapshot;
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (first) {
      first = false;
      continue;
    }
    const auto comma = line.find(',');
    if (comma == std::string::npos || comma == 0 || comma + 1 >= line.size()) {
      continue;
    }

    const std::string instance = line.substr(0, comma);
    const std::string value = line.substr(comma + 1);
    try {
      snapshot[instance] = std::stod(value);
    } catch (...) {
      // Ignore malformed values; test will fail on missing entries.
    }
  }
  return snapshot;
}

}  // namespace

TEST_CASE("Heuristic stage UB gate: ILS snapshot matches on SPPRCLIB + Roberti",
          "[heuristic][ub_gate][stage][slow]") {
  constexpr int kRestarts = 16;

  const auto instances = collect_heuristic_gate_instances();
  REQUIRE_FALSE(instances.empty());

  const auto expected = load_expected_snapshot();
  REQUIRE_FALSE(expected.empty());

  int matched = 0;
  std::set<std::string> seen;

  std::cout << "\nHeuristic stage UB snapshot gate (restarts=" << kRestarts
            << ")\n";
  std::cout << std::left << std::setw(32) << "Instance" << std::right
            << std::setw(14) << "Expected" << std::setw(14) << "Observed"
            << std::setw(14) << "Delta" << std::setw(14) << "Tol\n";

  for (const auto& path : instances) {
    const std::string key = path.filename().string();
    seen.insert(key);

    const auto it = expected.find(key);
    INFO("instance=" << key);
    REQUIRE(it != expected.end());

    auto prob = cptp::io::load(path);
    auto ils = cptp::heuristic::build_warm_start(prob, kRestarts);
    REQUIRE(std::isfinite(ils.objective));

    const double expected_obj = it->second;
    const double delta = ils.objective - expected_obj;
    // Snapshot values are stored from compact textual output (~6 significant
    // digits).
    const double tol = 1e-3 + 1e-5 * std::abs(expected_obj);

    std::cout << std::left << std::setw(32) << key << std::right
              << std::setw(14) << expected_obj << std::setw(14) << ils.objective
              << std::setw(14) << delta << std::setw(14) << tol << "\n";

    REQUIRE(std::abs(delta) <= tol);
    matched++;
  }

  for (const auto& [name, _] : expected) {
    INFO("snapshot_entry=" << name);
    REQUIRE(seen.contains(name));
  }

  std::cout << "Summary: matched=" << matched
            << ", instances=" << instances.size()
            << ", snapshot_rows=" << expected.size() << "\n";
  REQUIRE(matched == static_cast<int>(instances.size()));
  REQUIRE(expected.size() == instances.size());
}
