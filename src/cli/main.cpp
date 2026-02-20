#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "core/io.h"
#include "model/model.h"

static void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <instance> [--<highs_option> <value> ...]\n"
              << "\nAll options are forwarded to HiGHS. Common ones:\n"
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

    // Everything after the instance path is --key value pairs for HiGHS
    cptp::SolverOptions options;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.starts_with("--") && i + 1 < argc) {
            options.emplace_back(arg.substr(2), argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    try {
        auto problem = cptp::io::load(instance_path);
        std::cout << "Instance: " << problem.name
                  << " (" << problem.num_nodes() << " nodes, "
                  << problem.num_edges() << " edges)\n";

        cptp::Model model;
        model.set_problem(std::move(problem));

        auto result = model.solve(options);

        // Print solution
        std::cout << "\nSolution: ";
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

        if (result.total_cuts > 0) {
            std::cout << "User cuts: " << result.total_cuts
                      << " (" << result.separation_rounds << " rounds)\n";
            std::cout << std::fixed << std::setprecision(3);
            for (const auto& [name, stats] : result.separator_stats) {
                std::cout << "  " << std::setw(10) << std::left << name
                          << std::right
                          << std::setw(6) << stats.cuts_added << " cuts"
                          << std::setw(6) << stats.rounds_called << " rounds"
                          << std::setw(8) << stats.time_seconds << "s\n";
            }
        }

        return result.has_solution() ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
