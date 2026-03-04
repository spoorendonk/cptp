#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "core/io.h"
#include "model/model.h"

static void print_usage(const char* prog) {
    std::cerr << "cptp-solve — Capacitated Profitable Tour Problem solver\n\n"
              << "Usage: " << prog
              << " <instance> [--source <node>] [--target <node>] [--<highs_option> <value> ...]\n"
              << "\nOptions:\n"
              << "  --source <node>        Source node (overrides file value)\n"
              << "  --target <node>        Target node (overrides file value)\n"
              << "                         If source != target, solves s-t path instead of tour\n"
              << "  --cutoff <value>       Known UB for prove-only mode (skip finding incumbent)\n"
              << "  --disable_heuristics true  Disable warm-start and HiGHS MIP heuristics\n"
              << "\nSolver options:\n"
              << "  --rc_fixing <strategy> Lagrangian reduced-cost fixing strategy:\n"
              << "                         off (default) | root_only | on_ub_improvement | periodic | adaptive\n"
              << "  --rc_fixing_interval N Interval for periodic strategy (default 100)\n"
              << "  --rc_fixing_to_one true  Enable fix-to-1 for node variables (expensive)\n"
              << "  --deterministic_work_units N  Global cap for non-HiGHS heuristic work\n"
              << "                         (0 = uncapped/count-only, default)\n"
              << "  --heuristic_deterministic_restarts N  Fixed restarts per callback in deterministic mode\n"
              << "  --max_concurrent_solves N  Limit number of concurrent model solves in-process\n"
              << "                             (0 = no per-call cap; default)\n"
              << "  --heuristic_node_interval N  Heuristic callback cadence in B&B nodes (default 200)\n"
              << "  --enable_sec/--enable_rci/--enable_multistar/--enable_comb/--enable_spi true|false\n"
              << "                         Enable or disable cut families\n"
              << "  --max_cuts_<family> N  Per-family cut cap (sec,rci,multistar,comb,rglm,spi)\n"
              << "  --min_violation_<family> X  Per-family min violation filter\n"
              << "  --edge_elimination true/false  Enable preprocessing edge elimination (default true)\n"
              << "  --edge_elimination_nodes true/false  Also fix nodes if all incident edges are eliminated\n"
              << "  --preproc_stage1_bounds <mode>  Stage-1 bounds backend: two_cycle | ng_dssr | auto (default two_cycle)\n"
              << "  --workflow_dump true/false  Print startup/solve DAG wiring (default false)\n"
              << "  --branch_hyper_sb_max_depth/iter_limit/min_reliable/max_candidates <int>\n"
              << "                         Hyperplane strong-branching tuning\n"
              << "\nAll other options are forwarded to HiGHS. Common ones:\n"
              << "  --time_limit <sec>     Time limit (default 600)\n"
              << "  --threads <n>          Thread count\n"
              << "  --output_flag false    Suppress HiGHS output\n"
              << "\nSee HiGHS documentation for the full list.\n";
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

    std::filesystem::path instance_path = first_arg;

    // Parse arguments: extract --source/--target, forward rest to HiGHS
    rcspp::SolverOptions options;
    int override_source = -1;
    int override_target = -1;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--source" && i + 1 < argc) {
            override_source = std::stoi(argv[++i]);
        } else if (arg == "--target" && i + 1 < argc) {
            override_target = std::stoi(argv[++i]);
        } else if (arg.starts_with("--") && i + 1 < argc) {
            options.emplace_back(arg.substr(2), argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    try {
        auto problem = rcspp::io::load(instance_path);

        // Apply CLI overrides for source/target
        if (override_source >= 0 || override_target >= 0) {
            int32_t src = (override_source >= 0) ? override_source : problem.source();
            int32_t tgt = (override_target >= 0) ? override_target : problem.target();

            // Rebuild problem with new source/target
            problem.build(problem.num_nodes(),
                          // We need the edges — rebuild via the graph
                          // Actually, set_problem on Model handles this already.
                          // For CLI override, we re-build the problem.
                          [&] {
                              std::vector<rcspp::Edge> edges;
                              const auto& g = problem.graph();
                              for (auto e : g.edges())
                                  edges.push_back({g.edge_source(e), g.edge_target(e)});
                              return edges;
                          }(),
                          problem.edge_costs(),
                          problem.profits(),
                          problem.demands(),
                          problem.capacity(),
                          src, tgt);
        }

        std::cout << "Instance: " << problem.name
                  << " (" << problem.num_nodes() << " nodes, "
                  << problem.num_edges() << " edges";
        if (problem.is_tour()) {
            std::cout << ", tour from depot " << problem.source();
        } else {
            std::cout << ", path " << problem.source() << " -> " << problem.target();
        }
        std::cout << ")\n";

        rcspp::Model model;
        model.set_problem(std::move(problem));

        auto result = model.solve(options);

        // Print solution
        std::cout << "\n" << (model.problem().is_tour() ? "Tour" : "Path") << ": ";
        for (size_t i = 0; i < result.tour.size(); ++i) {
            if (i > 0) std::cout << " -> ";
            std::cout << result.tour[i];
        }
        std::cout << "\n";

        std::cout << "Objective: " << result.objective
                  << "  Bound: " << result.bound
                  << "  Gap: " << (result.gap * 100.0) << "%"
                  << "  Time: " << result.time_seconds << "s"
                  << "  Nodes: " << result.nodes << "\n";

        if (!result.separator_stats.empty()) {
            std::cout << "User cuts: " << result.total_cuts
                      << " (" << result.separation_rounds << " rounds)\n";
            std::cout << std::fixed << std::setprecision(3);
            const std::array<std::string, 6> preferred = {
                "SEC", "RCI", "Multistar", "Comb", "RGLM", "SPI"
            };
            auto print_sep_row = [](const std::string& name,
                                    const rcspp::SeparatorStats& stats) {
                std::cout << "  " << std::setw(10) << std::left << name
                          << std::right
                          << std::setw(6) << stats.cuts_added << " cuts"
                          << std::setw(6) << stats.rounds_called << " rounds"
                          << std::setw(8) << stats.time_seconds << "s\n";
            };
            for (const auto& name : preferred) {
                auto it = result.separator_stats.find(name);
                if (it != result.separator_stats.end()) {
                    print_sep_row(it->first, it->second);
                }
            }
            for (const auto& [name, stats] : result.separator_stats) {
                if (std::find(preferred.begin(), preferred.end(), name) == preferred.end()) {
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
