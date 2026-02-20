#include "model/model.h"
#include "model/highs_bridge.h"

#include <iostream>
#include <thread>

#include "heuristic/warm_start.h"
#include "sep/sec_separator.h"
#include "sep/rci_separator.h"
#include "sep/multistar_separator.h"
#include "sep/comb_separator.h"

#include "util/timer.h"

namespace cptp {

Model::Model() = default;

void Model::set_problem(Problem prob) {
    problem_ = std::move(prob);
    built_ = true;
}

void Model::set_graph(int32_t num_nodes,
                      std::span<const Edge> edges,
                      std::span<const double> edge_costs) {
    num_nodes_ = num_nodes;
    edges_.assign(edges.begin(), edges.end());
    edge_costs_.assign(edge_costs.begin(), edge_costs.end());
    built_ = false;
}

void Model::set_depot(int32_t depot) {
    depot_ = depot;
    built_ = false;
}

void Model::set_profits(std::span<const double> profits) {
    profits_.assign(profits.begin(), profits.end());
    built_ = false;
}

void Model::add_capacity_resource(std::span<const double> demands, double limit) {
    Resource res;
    res.type = Resource::Type::Capacity;
    res.node_consumption.assign(demands.begin(), demands.end());
    res.limit = limit;
    resources_.push_back(std::move(res));
    built_ = false;
}

void Model::build_problem() {
    if (built_) return;

    // Use demands/capacity from first capacity resource, or defaults
    std::vector<double> demands(num_nodes_, 0.0);
    double capacity = 1e18;

    for (const auto& res : resources_) {
        if (res.type == Resource::Type::Capacity) {
            demands = res.node_consumption;
            capacity = res.limit;
            break;
        }
    }

    if (profits_.empty()) {
        profits_.assign(num_nodes_, 0.0);
    }

    problem_.build(num_nodes_, edges_, edge_costs_, profits_, demands,
                   capacity, depot_);
    built_ = true;
}

SolveResult Model::solve(const SolverOptions& options) {
    build_problem();

    Timer timer;

    Highs highs;
    // Our defaults (user options can override)
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads",
        static_cast<int>(std::thread::hardware_concurrency()));

    // Forward all user options to HiGHS
    for (const auto& [key, value] : options) {
        auto status = highs.setOptionValue(key, value);
        if (status != HighsStatus::kOk) {
            std::cerr << "Warning: HiGHS rejected option " << key
                      << " = " << value << "\n";
        }
    }

    HiGHSBridge bridge(problem_, highs);

    // Add separators
    bridge.add_separator(std::make_unique<sep::SECSeparator>());
    bridge.add_separator(std::make_unique<sep::RCISeparator>());
    bridge.add_separator(std::make_unique<sep::MultistarSeparator>());
    bridge.add_separator(std::make_unique<sep::CombSeparator>());

    bridge.build_formulation();
    bridge.install_separators();

    // Warm-start with greedy insertion heuristic
    // Scale budget: ~1s for 50+ nodes, less for tiny instances
    {
        Timer ws_timer;
        double budget_ms = std::min(500.0, std::max(10.0,
            static_cast<double>(problem_.num_nodes()) * 10.0));
        auto warm = heuristic::build_warm_start(problem_, budget_ms);
        HighsSolution start;
        start.value_valid = true;
        start.col_value = std::move(warm);
        highs.setSolution(start);
        std::cerr << "Warm-start heuristic: " << ws_timer.elapsed_seconds() << "s\n";
    }

    highs.run();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();

    return result;
}

}  // namespace cptp
