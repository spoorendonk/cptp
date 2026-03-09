#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "Highs.h"
#include "core/io.h"
#include "model/model.h"

static void print_usage(const char* prog) {
  // All parameter lines use a fixed 42-char left column so descriptions align.
  std::cerr
      << "cptp-solve — Capacitated Profitable Tour Problem solver\n\n"
      << "Usage: " << prog
      << " <instance> [--source <node>] [--target <node>] [--<option> <value> "
         "...]\n"
      //
      << "\nGeneral:\n"
      << "  --source <node>                          Source node (overrides "
         "file value)\n"
      << "  --target <node>                          Target node (overrides "
         "file value)\n"
      << "                                           If source != target, "
         "solves s-t path\n"
      << "  --cutoff <value>                         Known UB for prove-only "
         "mode\n"
      << "  --time_limit <sec>                       Wall-clock time limit "
         "(adjusted for preprocessing)\n"
      << "  --threads N                              HiGHS thread count\n"
      << "  --random_seed N                          Random seed for HiGHS and "
         "heuristics (default 0)\n"
      << "  --output_flag true/false                 HiGHS console output "
         "(default true)\n"
      //
      << "\nBounds & Preprocessing:\n"
      << "  --labeling_elim_bounds true/false        Capacity-aware labeling "
         "bounds (default true)\n"
      << "                                           When false, disables edge "
         "elim & propagation\n"
      << "  --labeling_max_queue_pops N              Max queue pops per "
         "labeling run\n"
      << "                                           (default 1000000, "
         "0=unlimited)\n"
      << "  --all_pairs_bounds true/false            All-pairs bounds for "
         "propagation + SPI\n"
      << "                                           (default false)\n"
      << "  --edge_elimination true/false            Edge elimination (default "
         "true, needs bounds)\n"
      << "  --edge_elimination_nodes true/false      Node fixing from "
         "eliminated edges\n"
      << "                                           (default true, needs "
         "bounds)\n"
      //
      << "\nPropagation:\n"
      << "  --bounds_propagation true/false          Bounds-based domain "
         "propagation\n"
      << "                                           (default true, needs "
         "bounds)\n"
      << "  --rc_fixing <strategy>                   Reduced-cost fixing "
         "strategy\n"
      << "                                           "
         "off|root_only|on_ub_improvement|periodic|adaptive\n"
      << "  --rc_fixing_interval N                   Interval for periodic "
         "strategy (default 100)\n"
      << "  --rc_fixing_to_one true/false            Fix-to-1 for node "
         "variables\n"
      << "                                           (default false, "
         "expensive)\n"
      //
      << "\nHeuristics:\n"
      << "  Warm-start (construction + local search):\n"
      << "  --heu_ws true/false                      Enable warm-start "
         "heuristic (default true)\n"
      << "  --heu_ws_ls_max_iter N                   Local search iteration "
         "limit (default 1000)\n"
      << "\n"
      << "  LP-guided (B&C callback):\n"
      << "  --heu_lpg true/false                     Enable LP-guided "
         "heuristic (default true)\n"
      << "  --heu_lpg_strategy N                     0=all, 1=LP-threshold, "
         "2=RINS, 3=neighborhood\n"
      << "  --heu_lpg_budget_ms N                    Time budget per callback "
         "in ms (default 20)\n"
      << "  --heu_lpg_node_interval N                Callback cadence in B&B "
         "nodes (default 200)\n"
      << "  --heu_lpg_deterministic_restarts N       Fixed restarts per "
         "callback\n"
      << "  --heu_lpg_edge_threshold X               LP-threshold edge cutoff "
         "(default 0.1)\n"
      << "  --heu_lpg_node_threshold X               LP-threshold node cutoff "
         "(default 0.5)\n"
      << "  --heu_lpg_lp_threshold X                 RINS LP edge cutoff "
         "(default 0.1)\n"
      << "  --heu_lpg_seed_threshold X               Neighborhood seed edge "
         "cutoff (default 0.3)\n"
      << "\n"
      << "  HiGHS sub-MIP:\n"
      << "  --heu_highs_submip_sec true/false        SEC in HiGHS sub-MIPs "
         "(default true)\n"
      //
      << "\nSeparation (global):\n"
      << "  --separation_interval N                  Separate every N-th "
         "callback (default 1)\n"
      << "  --separation_tol X                       Fractional separation "
         "tolerance (default 0.1)\n"
      << "  --max_cuts_per_round N                   Global cut cap: pool all "
         "cuts, sort by\n"
      << "                                           violation, add top N "
         "(default 10, 0=unlimited)\n"
      //
      << "\nSeparation (per cut family):\n"
      << "  --enable_<family> true/false             Enable/disable a cut "
         "family\n"
      << "\n"
      << "  --enable_sec true/false                  SEC (default true)\n"
      << "  --enable_rci true/false                  RCI (default true)\n"
      << "  --enable_multistar true/false            Multistar (default true)\n"
      << "  --enable_comb true/false                 Comb (default false)\n"
      << "  --enable_rglm true/false                 RGLM (default false)\n"
      << "\n"
      << "  SPI (requires --all_pairs_bounds true):\n"
      << "  --enable_spi true/false                  SPI (default false)\n"
      << "  --spi_max_dp_size N                      Largest node set for "
         "exact Held-Karp\n"
      << "                                           lower bound; O(2^N·N²) "
         "(default 15)\n"
      << "  --spi_max_seeds N                        Infeasible pairs to grow "
         "into larger\n"
      << "                                           cuts, most-violated first "
         "(default 50)\n"
      << "  --spi_max_grow_candidates N              High-y nodes tried when "
         "extending a\n"
      << "                                           seed pair (default 40)\n"
      //
      << "\nBranching:\n"
      << "  --branch_hyper <mode>                    Hyperplane branching "
         "mode\n"
      << "                                           "
         "off|pairs|clusters|demand|cardinality|all\n"
      << "  --branch_hyper_sb_max_depth N            Strong-branching max "
         "depth\n"
      << "  --branch_hyper_sb_iter_limit N           Strong-branching "
         "iteration limit\n"
      << "  --branch_hyper_sb_min_reliable N         Strong-branching min "
         "reliable\n"
      << "  --branch_hyper_sb_max_candidates N       Strong-branching max "
         "candidates\n"
      //
      << "\nHiGHS:\n"
      << "  Unrecognized --key value pairs are forwarded to HiGHS.\n"
      << "  Presolve is always disabled (cptp requires stable column "
         "mapping).\n"
      << "  --highs_help                             Print all HiGHS options "
         "and exit\n";
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    print_usage(argv[0]);
    return 1;
  }

  std::string first_arg = argv[1];
  if (first_arg == "-h" || first_arg == "--help") {
    print_usage(argv[0]);
    return 0;
  }
  if (first_arg == "--highs_help") {
    Highs highs;
    highs.writeOptions("");
    return 0;
  }

  std::filesystem::path instance_path = first_arg;

  // Parse arguments: extract --source/--target, forward rest to HiGHS
  cptp::SolverOptions options;
  int override_source = -1;
  int override_target = -1;

  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--source" && i + 1 < argc) {
      override_source = std::stoi(argv[++i]);
    } else if (arg == "--target" && i + 1 < argc) {
      override_target = std::stoi(argv[++i]);
    } else if (arg == "--highs_help") {
      Highs highs;
      highs.writeOptions("");
      return 0;
    } else if (arg.starts_with("--") && i + 1 < argc) {
      options.emplace_back(arg.substr(2), argv[++i]);
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  try {
    auto problem = cptp::io::load(instance_path);

    // Apply CLI overrides for source/target
    if (override_source >= 0 || override_target >= 0) {
      int32_t src = (override_source >= 0) ? override_source : problem.source();
      int32_t tgt = (override_target >= 0) ? override_target : problem.target();

      // Rebuild problem with new source/target
      problem.build(
          problem.num_nodes(),
          // We need the edges — rebuild via the graph
          // Actually, set_problem on Model handles this already.
          // For CLI override, we re-build the problem.
          [&] {
            std::vector<cptp::Edge> edges;
            const auto& g = problem.graph();
            for (auto e : g.edges())
              edges.push_back({g.edge_source(e), g.edge_target(e)});
            return edges;
          }(),
          problem.edge_costs(), problem.profits(), problem.demands(),
          problem.capacity(), src, tgt);
    }

    std::cout << "Instance: " << problem.name << " (" << problem.num_nodes()
              << " nodes, " << problem.num_edges() << " edges";
    if (problem.is_tour()) {
      std::cout << ", tour from depot " << problem.source();
    } else {
      std::cout << ", path " << problem.source() << " -> " << problem.target();
    }
    std::cout << ")\n";

    cptp::Model model;
    model.set_problem(std::move(problem));

    auto result = model.solve(options);

    // Print solution
    std::cout << "\n" << (model.problem().is_tour() ? "Tour" : "Path") << ": ";
    for (size_t i = 0; i < result.tour.size(); ++i) {
      if (i > 0) std::cout << " -> ";
      std::cout << result.tour[i];
    }
    std::cout << "\n";

    std::cout << std::setprecision(15) << "Objective: " << result.objective
              << "  Bound: " << result.bound
              << "  Gap: " << (result.gap * 100.0) << "%"
              << "  Time: " << result.time_seconds << "s"
              << "  Nodes: " << result.nodes << "\n";

    if (!result.separator_stats.empty()) {
      std::cout << "User cuts: " << result.total_cuts << " ("
                << result.separation_rounds << " rounds)\n";
      std::cout << std::fixed << std::setprecision(3);
      const std::array<std::string, 6> preferred = {"SEC",  "RCI",  "Multistar",
                                                    "Comb", "RGLM", "SPI"};
      auto print_sep_row = [](const std::string& name,
                              const cptp::SeparatorStats& stats) {
        std::cout << "  " << std::setw(10) << std::left << name << std::right
                  << std::setw(6) << stats.cuts_added << " cuts" << std::setw(6)
                  << stats.rounds_called << " rounds" << std::setw(8)
                  << stats.time_seconds << "s\n";
      };
      for (const auto& name : preferred) {
        auto it = result.separator_stats.find(name);
        if (it != result.separator_stats.end()) {
          print_sep_row(it->first, it->second);
        }
      }
      for (const auto& [name, stats] : result.separator_stats) {
        if (std::find(preferred.begin(), preferred.end(), name) ==
            preferred.end()) {
          print_sep_row(name, stats);
        }
      }
    }

    return result.has_solution() ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
