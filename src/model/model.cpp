#include "model/model.h"
#include "model/highs_bridge.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cmath>
#include <mutex>
#include <set>
#include <thread>

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>

#include "heuristic/primal_heuristic.h"
#include "model/async_incumbent.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/ng_labeling.h"
#include "preprocess/shared_bounds.h"
#include "mip/HighsUserSeparator.h"
#include "mip/HighsMipSolverData.h"
#include "sep/sec_separator.h"
#include "sep/rci_separator.h"
#include "sep/multistar_separator.h"
#include "sep/comb_separator.h"
#include "sep/rglm_separator.h"
#include "sep/spi_separator.h"

#include "util/timer.h"

namespace rcspp {

namespace {

class SolveConcurrencyGuard {
 public:
    explicit SolveConcurrencyGuard(int32_t limit)
        : active_(false),
          requested_limit_(limit > 0 ? limit : kNoCap),
          limited_(requested_limit_ != kNoCap) {
        std::unique_lock<std::mutex> lock(mutex_);
        // Register finite limits before waiting so pending strict requests
        // immediately constrain subsequent arrivals.
        if (limited_) {
            limit_it_ = active_limits_.insert(requested_limit_);
        }
        cv_.wait(lock, [&] {
            return active_solves_ < std::min(requested_limit_, current_cap_locked());
        });
        ++active_solves_;
        active_ = true;
    }

    ~SolveConcurrencyGuard() {
        if (!active_) return;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (limited_) {
                active_limits_.erase(limit_it_);
            }
            --active_solves_;
        }
        cv_.notify_all();
    }

 private:
    static constexpr int32_t kNoCap = std::numeric_limits<int32_t>::max();
    static int32_t current_cap_locked() {
        return active_limits_.empty() ? kNoCap : *active_limits_.begin();
    }

    bool active_;
    int32_t requested_limit_;
    bool limited_;
    std::multiset<int32_t>::iterator limit_it_;
    static std::mutex mutex_;
    static std::condition_variable cv_;
    static int32_t active_solves_;
    static std::multiset<int32_t> active_limits_;
};

std::mutex SolveConcurrencyGuard::mutex_;
std::condition_variable SolveConcurrencyGuard::cv_;
int32_t SolveConcurrencyGuard::active_solves_ = 0;
std::multiset<int32_t> SolveConcurrencyGuard::active_limits_;

heuristic::HeuristicResult build_ng_path_start(const Problem& problem,
                                               std::span<const int32_t> path) {
    heuristic::HeuristicResult result;
    result.objective = std::numeric_limits<double>::infinity();
    if (path.size() < 2) return result;

    const auto& graph = problem.graph();
    const int32_t m = problem.num_edges();
    const int32_t n = problem.num_nodes();

    std::vector<double> col_values(static_cast<size_t>(m + n), 0.0);

    for (int32_t node : path) {
        if (node < 0 || node >= n) return result;
        col_values[static_cast<size_t>(m + node)] = 1.0;
    }

    double objective = 0.0;
    for (size_t k = 0; k + 1 < path.size(); ++k) {
        const int32_t u = path[k];
        const int32_t v = path[k + 1];
        int32_t edge = -1;
        for (auto e : graph.incident_edges(u)) {
            if (graph.other_endpoint(e, u) == v) {
                edge = e;
                break;
            }
        }
        if (edge < 0) return result;
        col_values[edge] = 1.0;
        objective += problem.edge_cost(edge);
    }

    for (int32_t i = 0; i < n; ++i) {
        if (col_values[static_cast<size_t>(m + i)] > 0.5) {
            objective -= problem.profit(i);
        }
    }

    result.col_values = std::move(col_values);
    result.objective = objective;
    return result;
}

}  // namespace

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
    std::string branch_hyper = "off";
    bool output_flag = true;
    bool deterministic = true;
    int64_t heuristic_node_interval = 200;
    bool heuristic_async_injection = true;
    int32_t max_concurrent_solves = 0;
    bool enable_sec = true;
    bool enable_rci = true;
    bool enable_multistar = true;
    bool enable_comb = true;
    bool enable_spi = true;
    bool edge_elimination = true;
    bool edge_elimination_nodes = true;
    double cutoff = std::numeric_limits<double>::infinity();
    bool disable_heuristics = false;
    bool dssr_async = true;
    int32_t hyper_sb_max_depth = 0;
    int32_t hyper_sb_iter_limit = 100;
    int32_t hyper_sb_min_reliable = 4;
    int32_t hyper_sb_max_candidates = 3;
    preprocess::ng::DssrOptions dssr_options;
    RCFixingSettings rc_settings;
    int32_t max_cuts_sec = -1;
    int32_t max_cuts_rci = -1;
    int32_t max_cuts_multistar = -1;
    int32_t max_cuts_comb = -1;
    int32_t max_cuts_rglm = -1;
    int32_t max_cuts_spi = -1;
    double min_violation_sec = -1.0;
    double min_violation_rci = -1.0;
    double min_violation_multistar = -1.0;
    double min_violation_comb = -1.0;
    double min_violation_rglm = -1.0;
    double min_violation_spi = -1.0;
    for (const auto& [key, value] : options) {
        if (key == "presolve") {
            // Our callback/cut plumbing assumes original column space and
            // stable row/col mappings; keep presolve disabled for now.
            if (value != "off" && value != "false" && value != "0") {
                logger_.log("Warning: option presolve={} requested, but rcspp-bac currently forces presolve=off", value);
            }
            continue;
        }
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
        if (key == "heuristic_node_interval") {
            heuristic_node_interval = std::stoll(value);
            continue;
        }
        if (key == "heuristic_async_injection") {
            heuristic_async_injection = (value == "true" || value == "1");
            continue;
        }
        if (key == "deterministic") {
            deterministic = (value == "true" || value == "1");
            continue;
        }
        if (key == "branch_hyper") {
            branch_hyper = value;
            continue;
        }
        if (key == "cutoff") {
            cutoff = std::stod(value);
            continue;
        }
        if (key == "max_concurrent_solves") {
            max_concurrent_solves = std::stoi(value);
            continue;
        }
        if (key == "enable_sec") {
            enable_sec = (value == "true" || value == "1");
            continue;
        }
        if (key == "enable_rci") {
            enable_rci = (value == "true" || value == "1");
            continue;
        }
        if (key == "enable_multistar") {
            enable_multistar = (value == "true" || value == "1");
            continue;
        }
        if (key == "enable_comb") {
            enable_comb = (value == "true" || value == "1");
            continue;
        }
        if (key == "enable_spi") {
            enable_spi = (value == "true" || value == "1");
            continue;
        }
        if (key == "edge_elimination") {
            edge_elimination = (value == "true" || value == "1");
            continue;
        }
        if (key == "edge_elimination_nodes") {
            edge_elimination_nodes = (value == "true" || value == "1");
            continue;
        }
        if (key == "disable_heuristics") {
            disable_heuristics = (value == "true" || value == "1");
            continue;
        }
        if (key == "dssr_async") {
            dssr_async = (value == "true" || value == "1");
            continue;
        }
        if (key == "ng_initial_size") {
            dssr_options.initial_ng_size = std::stoi(value);
            continue;
        }
        if (key == "ng_max_size") {
            dssr_options.max_ng_size = std::stoi(value);
            continue;
        }
        if (key == "ng_dssr_iters") {
            dssr_options.dssr_iterations = std::stoi(value);
            continue;
        }
        if (key == "ng_label_budget") {
            dssr_options.max_labels_per_node = std::stoi(value);
            continue;
        }
        if (key == "ng_simd") {
            dssr_options.enable_simd = (value == "true" || value == "1");
            continue;
        }
        if (key == "rc_fixing") {
            if (value == "root_only") rc_settings.strategy = RCFixingStrategy::root_only;
            else if (value == "on_ub_improvement") rc_settings.strategy = RCFixingStrategy::on_ub_improvement;
            else if (value == "periodic") rc_settings.strategy = RCFixingStrategy::periodic;
            else if (value == "adaptive" || value == "auto") rc_settings.strategy = RCFixingStrategy::adaptive;
            else rc_settings.strategy = RCFixingStrategy::off;
            continue;
        }
        if (key == "rc_fixing_interval") {
            rc_settings.periodic_interval = std::stoi(value);
            continue;
        }
        if (key == "rc_fixing_to_one") {
            rc_settings.fix_to_one = (value == "true" || value == "1");
            continue;
        }
        if (key == "branch_hyper_sb_max_depth") {
            hyper_sb_max_depth = std::stoi(value);
            continue;
        }
        if (key == "branch_hyper_sb_iter_limit") {
            hyper_sb_iter_limit = std::stoi(value);
            continue;
        }
        if (key == "branch_hyper_sb_min_reliable") {
            hyper_sb_min_reliable = std::stoi(value);
            continue;
        }
        if (key == "branch_hyper_sb_max_candidates") {
            hyper_sb_max_candidates = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_sec") {
            max_cuts_sec = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_rci") {
            max_cuts_rci = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_multistar") {
            max_cuts_multistar = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_comb") {
            max_cuts_comb = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_rglm") {
            max_cuts_rglm = std::stoi(value);
            continue;
        }
        if (key == "max_cuts_spi") {
            max_cuts_spi = std::stoi(value);
            continue;
        }
        if (key == "min_violation_sec") {
            min_violation_sec = std::stod(value);
            continue;
        }
        if (key == "min_violation_rci") {
            min_violation_rci = std::stod(value);
            continue;
        }
        if (key == "min_violation_multistar") {
            min_violation_multistar = std::stod(value);
            continue;
        }
        if (key == "min_violation_comb") {
            min_violation_comb = std::stod(value);
            continue;
        }
        if (key == "min_violation_rglm") {
            min_violation_rglm = std::stod(value);
            continue;
        }
        if (key == "min_violation_spi") {
            min_violation_spi = std::stod(value);
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

    dssr_options.initial_ng_size = std::max<int32_t>(1, dssr_options.initial_ng_size);
    dssr_options.max_ng_size = std::max<int32_t>(dssr_options.initial_ng_size,
                                                 dssr_options.max_ng_size);
    dssr_options.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations);
    dssr_options.max_labels_per_node = std::max<int32_t>(1, dssr_options.max_labels_per_node);
    heuristic_node_interval = std::max<int64_t>(1, heuristic_node_interval);
    max_concurrent_solves = std::max<int32_t>(0, max_concurrent_solves);
    hyper_sb_iter_limit = std::max<int32_t>(1, hyper_sb_iter_limit);
    hyper_sb_min_reliable = std::max<int32_t>(1, hyper_sb_min_reliable);
    hyper_sb_max_candidates = std::max<int32_t>(1, hyper_sb_max_candidates);
    hyper_sb_max_depth = std::max<int32_t>(0, hyper_sb_max_depth);

    SolveConcurrencyGuard solve_guard(max_concurrent_solves);

    // When cutoff is provided, set HiGHS objective_bound for node pruning.
    // This sets the initial upper_limit in the MIP solver, enabling pruning
    // of nodes whose LP bound exceeds the cutoff.
    // When heuristics are disabled, turn off HiGHS internal MIP heuristics.
    if (cutoff < std::numeric_limits<double>::infinity()) {
        highs.setOptionValue("objective_bound", cutoff);
    }
    if (disable_heuristics) {
        highs.setOptionValue("mip_heuristic_effort", 0.0);
    }
    // Enforce disabled presolve even if user passed --presolve.
    highs.setOptionValue("presolve", "off");

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
    heuristic::HeuristicResult warm_start{{}, std::numeric_limits<double>::infinity()};
    double warm_start_ub = std::numeric_limits<double>::infinity();
    int32_t ng_size_used = 1;
    double ng_path_ub = std::numeric_limits<double>::infinity();
    std::vector<int32_t> ng_path_nodes;

    // Use cutoff as UB when provided (known optimal for prove-only mode)
    if (cutoff < std::numeric_limits<double>::infinity()) {
        warm_start_ub = cutoff;
    }

    {
        Timer preproc_timer;

        if (all_pairs_propagation) {
            // All-pairs labeling runs in a parallel slot alongside warm-start.
            int num_restarts = std::clamp(n, 20, 200);
            double budget_ms = deterministic ? 0.0
                : std::max(10.0, static_cast<double>(n) * 10.0);
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
            if (!disable_heuristics) {
                tg.run([&] {
                    warm_start = heuristic::build_initial_solution(
                        problem_, num_restarts, budget_ms);
                });
            }
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
            // Default: source/target DSSR-ng labeling with warm-start
            int num_restarts = std::clamp(n, 20, 200);
            double budget_ms = deterministic ? 0.0
                : std::min(500.0, std::max(10.0,
                    static_cast<double>(n) * 10.0));

            tbb::task_group tg;
            tg.run([&] {
                auto bounds = preprocess::ng::compute_bounds(
                    problem_, source, target, dssr_options);
                fwd_bounds = std::move(bounds.fwd);
                bwd_bounds = std::move(bounds.bwd);
                ng_size_used = bounds.ng_size;
                if (bounds.elementary_path_found) {
                    ng_path_ub = bounds.elementary_path_cost;
                    ng_path_nodes = std::move(bounds.elementary_path);
                }
            });
            if (!disable_heuristics) {
                tg.run([&] {
                    warm_start = heuristic::build_initial_solution(
                        problem_, num_restarts, budget_ms);
                });
            }
            tg.wait();
        }

        // Heuristic UB only used when no cutoff was provided
        if (cutoff >= std::numeric_limits<double>::infinity()) {
            warm_start_ub = warm_start.objective;
        }
        if (ng_path_ub < warm_start_ub) {
            warm_start_ub = ng_path_ub;
        }
        if (!ng_path_nodes.empty()) {
            auto ng_start = build_ng_path_start(problem_, ng_path_nodes);
            if (!ng_start.col_values.empty() &&
                ng_start.objective < warm_start.objective - 1e-9) {
                warm_start = std::move(ng_start);
                warm_start_ub = std::min(warm_start_ub, warm_start.objective);
            }
        }
        logger_.log("Preprocessing: {}s, UB={}{}", preproc_timer.elapsed_seconds(),
                    warm_start_ub,
                    (cutoff < std::numeric_limits<double>::infinity() ? " (cutoff)" : ""));
    }

    HiGHSBridge bridge(problem_, highs, logger_, separation_tol);
    bridge.set_separation_interval(separation_interval);
    bridge.set_max_cuts_per_separator(max_cuts_per_sep);
    bridge.set_submip_separation(submip_separation);
    bridge.set_upper_bound(warm_start_ub);
    bridge.set_rc_fixing(rc_settings);
    bridge.set_edge_elimination(edge_elimination);
    bridge.set_edge_elimination_nodes(edge_elimination_nodes);
    bridge.set_heuristic_node_interval(heuristic_node_interval);
    bridge.set_heuristic_async_injection(heuristic_async_injection);
    if (max_cuts_sec >= 0) bridge.set_separator_max_cuts("SEC", max_cuts_sec);
    if (max_cuts_rci >= 0) bridge.set_separator_max_cuts("RCI", max_cuts_rci);
    if (max_cuts_multistar >= 0) bridge.set_separator_max_cuts("Multistar", max_cuts_multistar);
    if (max_cuts_comb >= 0) bridge.set_separator_max_cuts("Comb", max_cuts_comb);
    if (max_cuts_rglm >= 0) bridge.set_separator_max_cuts("RGLM", max_cuts_rglm);
    if (max_cuts_spi >= 0) bridge.set_separator_max_cuts("SPI", max_cuts_spi);
    if (min_violation_sec >= 0.0) bridge.set_separator_min_violation("SEC", min_violation_sec);
    if (min_violation_rci >= 0.0) bridge.set_separator_min_violation("RCI", min_violation_rci);
    if (min_violation_multistar >= 0.0) bridge.set_separator_min_violation("Multistar", min_violation_multistar);
    if (min_violation_comb >= 0.0) bridge.set_separator_min_violation("Comb", min_violation_comb);
    if (min_violation_rglm >= 0.0) bridge.set_separator_min_violation("RGLM", min_violation_rglm);
    if (min_violation_spi >= 0.0) bridge.set_separator_min_violation("SPI", min_violation_spi);

    auto shared_bounds = std::make_shared<preprocess::SharedBoundsStore>(
        preprocess::BoundSnapshot{
            .fwd = fwd_bounds,
            .bwd = bwd_bounds,
            .correction = correction,
            .ng_size = ng_size_used,
            .version = 0,
        });
    auto async_upper_bound = std::make_shared<std::atomic<double>>(warm_start_ub);
    auto async_incumbent_store = std::make_shared<model::AsyncIncumbentStore>();
    if (!warm_start.col_values.empty()) {
        async_incumbent_store->publish_if_better(
            std::vector<double>(warm_start.col_values.begin(), warm_start.col_values.end()),
            warm_start.objective);
    }
    bridge.set_shared_bounds_store(shared_bounds);
    bridge.set_async_upper_bound(async_upper_bound);
    bridge.set_async_incumbent_store(async_incumbent_store);

    bridge.set_labeling_bounds(std::move(fwd_bounds), std::move(bwd_bounds), correction);
    if (all_pairs_propagation) {
        bridge.set_all_pairs_bounds(std::move(all_pairs));
    }

    // Add separators (same default families for tour and s-t path).
    if (enable_sec)
        bridge.add_separator(std::make_unique<sep::SECSeparator>());
    if (enable_rci)
        bridge.add_separator(std::make_unique<sep::RCISeparator>());
    if (enable_multistar)
        bridge.add_separator(std::make_unique<sep::MultistarSeparator>());
    if (enable_rglm)
        bridge.add_separator(std::make_unique<sep::RGLMSeparator>());
    if (enable_comb)
        bridge.add_separator(std::make_unique<sep::CombSeparator>());
    if (all_pairs_propagation && enable_spi)
        bridge.add_separator(std::make_unique<sep::SPISeparator>());

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
        HighsUserSeparator::StrongBranchConfig sb_cfg;
        sb_cfg.max_depth = hyper_sb_max_depth;
        sb_cfg.iter_limit = hyper_sb_iter_limit;
        sb_cfg.min_reliable = hyper_sb_min_reliable;
        sb_cfg.max_sb_candidates = hyper_sb_max_candidates;
        HighsUserSeparator::setStrongBranchConfig(sb_cfg);

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
    bridge.set_heuristic_callback(heuristic_callback);
    bridge.set_heuristic_budget_ms(heuristic_budget_ms);
    bridge.set_heuristic_strategy(heuristic_strategy);
    bridge.install_heuristic_callback();

    // Pass initial solution to HiGHS (skip if heuristics disabled)
    if (!disable_heuristics) {
        HighsSolution start;
        start.value_valid = true;
        start.col_value = std::move(warm_start.col_values);
        HighsInt num_cols = highs.getNumCol();
        while (static_cast<HighsInt>(start.col_value.size()) < num_cols)
            start.col_value.push_back(0.0);
        highs.setSolution(start);
    }

    std::atomic<bool> stop_async_dssr{false};
    std::thread async_dssr_worker;
    if (dssr_async &&
        !all_pairs_propagation &&
        dssr_options.max_ng_size > ng_size_used) {
        async_dssr_worker = std::thread([&]() {
            const int32_t start_ng = std::max(ng_size_used + 1, dssr_options.initial_ng_size + 1);
            int32_t updates = 0;
            for (int32_t ng_size = start_ng;
                 ng_size <= dssr_options.max_ng_size &&
                 !stop_async_dssr.load(std::memory_order_relaxed);
                 ++ng_size) {
                auto stage_opts = dssr_options;
                stage_opts.initial_ng_size = ng_size;
                stage_opts.max_ng_size = ng_size;
                stage_opts.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations / 2);
                auto bounds = preprocess::ng::compute_bounds(problem_, source, target, stage_opts);
                bool changed = shared_bounds->publish_tightening(bounds.fwd, bounds.bwd, bounds.ng_size);
                if (changed) updates++;
                if (bounds.elementary_path_found) {
                    double prev = async_upper_bound->load(std::memory_order_relaxed);
                    while (bounds.elementary_path_cost + 1e-9 < prev &&
                           !async_upper_bound->compare_exchange_weak(
                               prev, bounds.elementary_path_cost, std::memory_order_relaxed)) {
                    }
                    if (!bounds.elementary_path.empty()) {
                        auto start = build_ng_path_start(problem_, bounds.elementary_path);
                        if (!start.col_values.empty()) {
                            async_incumbent_store->publish_if_better(
                                std::move(start.col_values), start.objective);
                        }
                    }
                }
            }
            (void)updates;
        });
    }

    highs.run();

    stop_async_dssr.store(true, std::memory_order_relaxed);
    if (async_dssr_worker.joinable()) {
        async_dssr_worker.join();
    }

    // Clear static callbacks/state to avoid leaking between solves
    HighsUserSeparator::clearCallback();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();

    return result;
}

}  // namespace rcspp
