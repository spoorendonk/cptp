#include "model/model.h"
#include "model/highs_bridge.h"

#include <algorithm>
#include <cmath>
#include <thread>

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>

#include "heuristic/warm_start.h"
#include "preprocess/edge_elimination.h"
#include "mip/HighsUserSeparator.h"
#include "mip/HighsMipSolverData.h"
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
    std::string branch_hyper = "off";
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
        if (key == "branch_hyper") {
            branch_hyper = value;
            continue;
        }
        if (key == "output_flag") {
            output_flag = (value == "true" || value == "1");
            continue;
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

    // Run preprocessing and warm-start in parallel:
    // - forward_labeling (source = depot)
    // - backward_labeling (target = depot, same as forward for undirected)
    // - build_warm_start
    const int32_t source = problem_.source();
    const int32_t target = problem_.target();
    // correction = profit(source) when s == t (tour: depot profit double-subtracted)
    double correction = problem_.is_tour() ? problem_.profit(source) : 0.0;

    const int32_t n = problem_.num_nodes();
    std::vector<double> fwd_bounds, bwd_bounds;
    std::vector<double> all_pairs;
    heuristic::WarmStartResult warm_start;
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
                warm_start = heuristic::build_warm_start(problem_, budget_ms);
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
                warm_start = heuristic::build_warm_start(problem_, budget_ms);
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

    // Hyperplane branching: dynamic constraint branching
    static constexpr std::string_view valid_modes[] = {
        "off", "pairs", "clusters", "demand", "cardinality", "all"
    };
    if (std::find(std::begin(valid_modes), std::end(valid_modes), branch_hyper)
            == std::end(valid_modes)) {
        logger_.log("Warning: unknown branch_hyper mode '{}', using 'off'",
                    branch_hyper);
        branch_hyper = "off";
    }
    if (branch_hyper != "off") {
        const int32_t n = bridge.num_nodes();
        const int32_t y_off = bridge.y_offset();
        const auto& graph = problem_.graph();

        std::vector<double> demands(n, 0.0);
        for (int32_t i = 0; i < n; ++i)
            demands[i] = problem_.demand(i);

        // Ryan-Foster pairs: for each non-depot node, pair with nearest neighbor
        struct Pair { int32_t i, j; };
        std::vector<Pair> rf_pairs;
        if (branch_hyper == "pairs" || branch_hyper == "all") {
            std::vector<bool> paired(n, false);
            for (int32_t i = 0; i < n; ++i) {
                if (i == problem_.source()) continue;
                double best_cost = std::numeric_limits<double>::infinity();
                int32_t best_j = -1;
                for (auto e : graph.incident_edges(i)) {
                    int32_t u = graph.edge_source(e);
                    int32_t v = graph.edge_target(e);
                    int32_t nb = (u == i) ? v : u;
                    if (nb == problem_.source()) continue;
                    double c = problem_.edge_cost(e);
                    if (c < best_cost) {
                        best_cost = c;
                        best_j = nb;
                    }
                }
                if (best_j >= 0) {
                    // Deduplicate: only add (min, max) pairs
                    int32_t lo = std::min(i, best_j);
                    int32_t hi = std::max(i, best_j);
                    if (!paired[lo] && !paired[hi]) {
                        rf_pairs.push_back({lo, hi});
                        paired[lo] = true;
                        paired[hi] = true;
                    }
                }
            }
        }

        // Greedy nearest-neighbor clusters of size ~4
        std::vector<std::vector<int32_t>> clusters;
        if (branch_hyper == "clusters" || branch_hyper == "all") {
            constexpr int32_t CLUSTER_SIZE = 4;
            std::vector<bool> assigned(n, false);
            assigned[problem_.source()] = true;

            // Sort nodes by demand (descending) to seed clusters
            std::vector<int32_t> nodes_by_demand;
            for (int32_t i = 0; i < n; ++i) {
                if (i != problem_.source())
                    nodes_by_demand.push_back(i);
            }
            std::sort(nodes_by_demand.begin(), nodes_by_demand.end(),
                [&](int32_t a, int32_t b) {
                    return demands[a] > demands[b];
                });

            for (int32_t seed : nodes_by_demand) {
                if (assigned[seed]) continue;
                std::vector<int32_t> cluster = {seed};
                assigned[seed] = true;

                // Add nearest unassigned neighbors
                while (static_cast<int32_t>(cluster.size()) < CLUSTER_SIZE) {
                    double best_cost = std::numeric_limits<double>::infinity();
                    int32_t best_nb = -1;
                    for (int32_t ci : cluster) {
                        for (auto e : graph.incident_edges(ci)) {
                            int32_t u = graph.edge_source(e);
                            int32_t v = graph.edge_target(e);
                            int32_t nb = (u == ci) ? v : u;
                            if (assigned[nb]) continue;
                            double c = problem_.edge_cost(e);
                            if (c < best_cost) {
                                best_cost = c;
                                best_nb = nb;
                            }
                        }
                    }
                    if (best_nb < 0) break;
                    cluster.push_back(best_nb);
                    assigned[best_nb] = true;
                }
                if (cluster.size() >= 2) {
                    clusters.push_back(std::move(cluster));
                }
            }
        }

        HighsUserSeparator::setBranchingCallback(
            [y_off, n, source = problem_.source(),
             demands = std::move(demands), branch_hyper,
             rf_pairs = std::move(rf_pairs),
             clusters = std::move(clusters)](
                const HighsMipSolver& mipsolver)
                -> std::vector<HighsUserSeparator::HyperplaneCandidate> {
                const auto& sol = mipsolver.mipdata_->lp.getSolution().col_value;
                std::vector<HighsUserSeparator::HyperplaneCandidate> result;

                if (branch_hyper == "pairs" || branch_hyper == "all") {
                    for (auto& [pi, pj] : rf_pairs) {
                        double yi = sol[y_off + pi], yj = sol[y_off + pj];
                        // Skip if either variable is near-integer
                        if (yi < 0.05 || yi > 0.95 || yj < 0.05 || yj > 0.95)
                            continue;
                        HighsUserSeparator::HyperplaneCandidate c;
                        c.indices = {static_cast<HighsInt>(y_off + pi),
                                     static_cast<HighsInt>(y_off + pj)};
                        c.values = {1.0, 1.0};
                        c.name = "pair_" + std::to_string(pi) + "_" + std::to_string(pj);
                        result.push_back(std::move(c));
                    }
                }
                if (branch_hyper == "clusters" || branch_hyper == "all") {
                    for (size_t k = 0; k < clusters.size(); ++k) {
                        auto& cluster = clusters[k];
                        // Cluster demand hyperplane: sum(demand_i * y_i) for i in cluster
                        HighsUserSeparator::HyperplaneCandidate c;
                        c.name = "cluster_" + std::to_string(k);
                        double total_demand = 0.0;
                        for (int32_t i : cluster) {
                            if (demands[i] > 0) {
                                c.indices.push_back(y_off + i);
                                c.values.push_back(demands[i]);
                                total_demand += demands[i] * sol[y_off + i];
                            }
                        }
                        if (!c.indices.empty()) {
                            double frac = total_demand - std::floor(total_demand);
                            if (frac > 0.01 && frac < 0.99)
                                result.push_back(std::move(c));
                        }
                    }
                }
                if (branch_hyper == "demand" || branch_hyper == "all") {
                    HighsUserSeparator::HyperplaneCandidate c;
                    c.name = "demand";
                    for (int32_t i = 0; i < n; ++i) {
                        if (i == source) continue;  // skip depot
                        if (demands[i] > 0) {
                            c.indices.push_back(y_off + i);
                            c.values.push_back(demands[i]);
                        }
                    }
                    if (!c.indices.empty()) result.push_back(std::move(c));
                }
                if (branch_hyper == "cardinality" || branch_hyper == "all") {
                    HighsUserSeparator::HyperplaneCandidate c;
                    c.name = "cardinality";
                    for (int32_t i = 0; i < n; ++i) {
                        if (i == source) continue;  // skip depot
                        c.indices.push_back(y_off + i);
                        c.values.push_back(1.0);
                    }
                    result.push_back(std::move(c));
                }
                return result;
            });
    }

    bridge.install_separators();
    bridge.install_propagator();

    // Pass warm-start solution to HiGHS
    {
        HighsSolution start;
        start.value_valid = true;
        start.col_value = std::move(warm_start.col_values);
        HighsInt num_cols = highs.getNumCol();
        while (static_cast<HighsInt>(start.col_value.size()) < num_cols)
            start.col_value.push_back(0.0);
        highs.setSolution(start);
    }

    highs.run();

    // Clear static callbacks/state to avoid leaking between solves
    HighsUserSeparator::clearCallback();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();

    return result;
}

}  // namespace rcspp
