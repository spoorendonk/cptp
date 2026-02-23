#include "model/model.h"
#include "model/highs_bridge.h"

#include <thread>

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>

#include "heuristic/primal_heuristic.h"
#include "preprocess/edge_elimination.h"
#include "sep/sec_separator.h"
#include "sep/rci_separator.h"
#include "sep/multistar_separator.h"
#include "sep/comb_separator.h"
#include "sep/rglm_separator.h"

#include "util/timer.h"

namespace rcspp {

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

void Model::set_source(int32_t source) {
    source_ = source;
    built_ = false;
}

void Model::set_target(int32_t target) {
    target_ = target;
    built_ = false;
}

void Model::set_depot(int32_t depot) {
    source_ = depot;
    target_ = depot;
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

    problem_.build(num_nodes_, edges_, std::move(edge_costs_),
                   std::move(profits_), std::move(demands),
                   capacity, source_, target_);
    built_ = true;
}

SolveResult Model::solve(const SolverOptions& options) {
    build_problem();

    logger_.log("rcspp-bac v0.1");

    Timer timer;

    Highs highs;
    // Our defaults (user options can override)
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads",
        static_cast<int>(std::thread::hardware_concurrency()));
    // Disable strong branching: trust pseudocosts immediately.
    // Benchmarks show pscost=0 or 2 outperforms default=8 on hard instances.
    highs.setOptionValue("mip_pscost_minreliable", 0);

    // Intercept our custom options, forward the rest to HiGHS
    int32_t separation_interval = 1;
    int32_t max_cuts_per_sep = 3;  // top-k most-violated per separator per round
    double separation_tol = sep::kDefaultFracTol;
    bool all_pairs_propagation = false;
    bool enable_rglm = false;
    bool submip_separation = true;
    bool heuristic_callback = true;
    double heuristic_budget_ms = 20.0;
    int heuristic_strategy = 0;
    bool output_flag = true;
    for (const auto& [key, value] : options) {
        if (key == "separation_interval") {
            separation_interval = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_per_separator") {
            max_cuts_per_sep = std::stoi(value);
            continue;
        }
        if (key == "separation_tol") {
            separation_tol = std::stod(value);
            continue;
        }
        if (key == "all_pairs_propagation") {
            all_pairs_propagation = (value == "true" || value == "1");
            continue;
        }
        if (key == "enable_rglm") {
            enable_rglm = (value == "true" || value == "1");
            continue;
        }
        if (key == "submip_separation") {
            submip_separation = (value == "true" || value == "1");
            continue;
        }
        if (key == "heuristic_callback") {
            heuristic_callback = (value == "true" || value == "1");
            continue;
        }
        if (key == "heuristic_budget_ms") {
            heuristic_budget_ms = std::stod(value);
            continue;
        }
        if (key == "heuristic_strategy") {
            heuristic_strategy = std::stoi(value);
            continue;
        }
        if (key == "output_flag") {
            output_flag = (value == "true" || value == "1");
        }
        auto status = highs.setOptionValue(key, value);
        if (status != HighsStatus::kOk) {
            logger_.log("Warning: HiGHS rejected option {} = {}", key, value);
        }
    }

    // Suppress HiGHS console output; forward through our logger if enabled
    highs.setOptionValue("output_flag", false);
    if (output_flag) {
        highs.setLogCallback(
            [](HighsLogType /*type*/, const char* message, void* data) {
                static_cast<Logger*>(data)->log(std::string_view{message});
            },
            &logger_);
    }

    // Run preprocessing and initial heuristic in parallel:
    // - forward_labeling (source = depot)
    // - backward_labeling (target = depot, same as forward for undirected)
    // - build_initial_solution
    const int32_t source = problem_.source();
    const int32_t target = problem_.target();
    // correction = profit(source) when s == t (tour: depot profit double-subtracted)
    double correction = problem_.is_tour() ? problem_.profit(source) : 0.0;

    const int32_t n = problem_.num_nodes();
    std::vector<double> fwd_bounds, bwd_bounds;
    std::vector<double> all_pairs;
    heuristic::HeuristicResult warm_start;
    double warm_start_ub = std::numeric_limits<double>::infinity();

    {
        Timer preproc_timer;

        if (all_pairs_propagation) {
            // All-pairs labeling: uncapped warm-start budget so it fills
            // the parallel slot (all-pairs typically takes longer).
            double budget_ms = std::max(10.0, static_cast<double>(n) * 10.0);
            constexpr double inf = std::numeric_limits<double>::infinity();
            all_pairs.assign(static_cast<size_t>(n) * n, inf);

            tbb::task_group tg;
            tg.run([&] {
                tbb::parallel_for(0, n, [&](int32_t s) {
                    auto row = preprocess::forward_labeling(problem_, s);
                    std::copy(row.begin(), row.end(),
                              all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
                });
            });
            tg.run([&] {
                warm_start = heuristic::build_initial_solution(problem_, budget_ms);
            });
            tg.wait();

            // Extract source row for fwd/bwd bounds
            fwd_bounds.assign(
                all_pairs.begin() + static_cast<ptrdiff_t>(source) * n,
                all_pairs.begin() + static_cast<ptrdiff_t>(source) * n + n);
            if (source == target) {
                bwd_bounds = fwd_bounds;
            } else {
                bwd_bounds.assign(
                    all_pairs.begin() + static_cast<ptrdiff_t>(target) * n,
                    all_pairs.begin() + static_cast<ptrdiff_t>(target) * n + n);
            }
        } else {
            // Default: depot-only labeling with capped warm-start budget
            double budget_ms = std::min(500.0, std::max(10.0,
                static_cast<double>(n) * 10.0));

            tbb::task_group tg;
            tg.run([&] {
                fwd_bounds = preprocess::forward_labeling(problem_, source);
            });
            if (source != target) {
                tg.run([&] {
                    bwd_bounds = preprocess::backward_labeling(problem_, target);
                });
            }
            tg.run([&] {
                warm_start = heuristic::build_initial_solution(problem_, budget_ms);
            });
            tg.wait();

            if (source == target) {
                bwd_bounds = fwd_bounds;
            }
        }

        warm_start_ub = warm_start.objective;
        logger_.log("Preprocessing: {}s, UB={}", preproc_timer.elapsed_seconds(),
                    warm_start_ub);
    }

    HiGHSBridge bridge(problem_, highs, logger_, separation_tol);
    bridge.set_separation_interval(separation_interval);
    bridge.set_max_cuts_per_separator(max_cuts_per_sep);
    bridge.set_submip_separation(submip_separation);
    bridge.set_upper_bound(warm_start_ub);
    bridge.set_labeling_bounds(std::move(fwd_bounds), std::move(bwd_bounds), correction);
    if (all_pairs_propagation) {
        bridge.set_all_pairs_bounds(std::move(all_pairs));
    }

    // Add separators
    bridge.add_separator(std::make_unique<sep::SECSeparator>());
    bridge.add_separator(std::make_unique<sep::RCISeparator>());
    bridge.add_separator(std::make_unique<sep::MultistarSeparator>());
    if (enable_rglm)
        bridge.add_separator(std::make_unique<sep::RGLMSeparator>());
    bridge.add_separator(std::make_unique<sep::CombSeparator>());

    bridge.build_formulation();
    bridge.install_separators();
    bridge.install_propagator();
    bridge.set_heuristic_callback(heuristic_callback);
    bridge.set_heuristic_budget_ms(heuristic_budget_ms);
    bridge.set_heuristic_strategy(heuristic_strategy);
    bridge.install_heuristic_callback();

    // Pass initial solution to HiGHS
    {
        HighsSolution start;
        start.value_valid = true;
        start.col_value = std::move(warm_start.col_values);
        highs.setSolution(start);
    }

    highs.run();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();

    return result;
}

}  // namespace rcspp
