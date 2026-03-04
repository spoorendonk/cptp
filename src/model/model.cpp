#include "model/model.h"
#include "model/highs_bridge.h"
#include "model/solve_concurrency_guard.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string_view>
#include <thread>

#ifdef __linux__
#include <unistd.h>
#endif

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include "heuristic/primal_heuristic.h"
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
        col_values[edge] += 1.0;
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

int32_t count_finite_bounds(std::span<const double> bounds) {
    return static_cast<int32_t>(std::count_if(
        bounds.begin(), bounds.end(),
        [](double v) { return std::isfinite(v); }));
}

bool has_feasible_heuristic(const heuristic::HeuristicResult& r) {
    return !r.col_values.empty() &&
           std::isfinite(r.objective) &&
           r.objective < std::numeric_limits<double>::max() / 2.0;
}

unsigned detect_physical_cores() {
#ifdef __linux__
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n > 0) return static_cast<unsigned>(n);
#endif
    unsigned logical = std::thread::hardware_concurrency();
    return (logical > 0) ? logical : 1u;
}

std::string join_tokens(const std::vector<std::string>& tokens,
                        std::string_view delim = " ") {
    if (tokens.empty()) return {};
    std::ostringstream oss;
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (i) oss << delim;
        oss << tokens[i];
    }
    return oss.str();
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

    logger_.log("cptp v0.1");

    Timer timer;

    Highs highs;
    // Our defaults (user options can override)
    highs.setOptionValue("presolve", "off");
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
    int64_t deterministic_work_units = 0;
    int32_t heuristic_deterministic_restarts = 32;
    int64_t heuristic_node_interval = 200;
    bool workflow_dump = false;
    bool warned_legacy_internal_option = false;
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
    std::string preproc_stage1_bounds = "two_cycle";
    int32_t preproc_fast_restarts = 12;
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
    std::map<std::string, std::string> accepted_highs_options;
    bool threads_option_explicit = false;
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
        if (key == "workflow_dump") {
            workflow_dump = (value == "true" || value == "1");
            continue;
        }
        if (key.rfind("_internal_", 0) == 0) {
            if (!warned_legacy_internal_option) {
                logger_.log("Warning: internal options are no longer supported and will be ignored");
                warned_legacy_internal_option = true;
            }
            continue;
        }
        if (key == "deterministic_work_units") {
            deterministic_work_units = std::stoll(value);
            continue;
        }
        if (key == "heuristic_deterministic_restarts") {
            heuristic_deterministic_restarts = std::stoi(value);
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
        if (key == "preproc_stage1_bounds") {
            std::string mode = value;
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (mode == "two_cycle" || mode == "ng_dssr" || mode == "auto") {
                preproc_stage1_bounds = mode;
            } else if (mode == "ng1") {
                preproc_stage1_bounds = "ng_dssr";
                logger_.log("Warning: preproc_stage1_bounds=ng1 is deprecated; use ng_dssr");
            } else {
                logger_.log("Warning: unknown preproc_stage1_bounds='{}' (expected two_cycle|ng_dssr|auto), using two_cycle",
                            value);
                preproc_stage1_bounds = "two_cycle";
            }
            continue;
        }
        if (key == "preproc_fast_restarts") {
            preproc_fast_restarts = std::stoi(value);
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
        if (key == "ng_simd") {
            dssr_options.enable_simd = (value == "true" || value == "1");
            continue;
        }
        if (key == "rc_fixing") {
            if (value == "root_only") rc_settings.strategy = RCFixingStrategy::root_only;
            else if (value == "on_ub_improvement") rc_settings.strategy = RCFixingStrategy::on_ub_improvement;
            else if (value == "periodic") rc_settings.strategy = RCFixingStrategy::periodic;
            else if (value == "adaptive") rc_settings.strategy = RCFixingStrategy::adaptive;
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
        } else {
            accepted_highs_options[key] = value;
            if (key == "threads") threads_option_explicit = true;
        }
    }

    dssr_options.initial_ng_size = std::max<int32_t>(0, dssr_options.initial_ng_size);
    dssr_options.max_ng_size = std::max<int32_t>(dssr_options.initial_ng_size,
                                                 dssr_options.max_ng_size);
    dssr_options.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations);
    preproc_fast_restarts = std::max<int32_t>(1, preproc_fast_restarts);
    heuristic_node_interval = std::max<int64_t>(1, heuristic_node_interval);
    deterministic_work_units = std::max<int64_t>(0, deterministic_work_units);
    heuristic_deterministic_restarts =
        std::max<int32_t>(1, heuristic_deterministic_restarts);
    max_concurrent_solves = std::max<int32_t>(0, max_concurrent_solves);
    hyper_sb_iter_limit = std::max<int32_t>(1, hyper_sb_iter_limit);
    hyper_sb_min_reliable = std::max<int32_t>(1, hyper_sb_min_reliable);
    hyper_sb_max_candidates = std::max<int32_t>(1, hyper_sb_max_candidates);
    hyper_sb_max_depth = std::max<int32_t>(0, hyper_sb_max_depth);
    HighsInt highs_threads_option = 0;
    if (highs.getOptionValue("threads", highs_threads_option) !=
        HighsStatus::kOk) {
        highs_threads_option = 0;
    }
    auto resolve_highs_auto_threads = []() -> int32_t {
        const unsigned hw = std::thread::hardware_concurrency();
        const unsigned resolved = (std::max(1u, hw) + 1u) / 2u;
        return static_cast<int32_t>(resolved);
    };
    int32_t preproc_thread_limit = 1;
    if (highs_threads_option > 0) {
        if (highs_threads_option >
            static_cast<HighsInt>(std::numeric_limits<int32_t>::max())) {
            preproc_thread_limit = std::numeric_limits<int32_t>::max();
        } else {
            preproc_thread_limit = static_cast<int32_t>(highs_threads_option);
        }
    } else {
        preproc_thread_limit = resolve_highs_auto_threads();
        // Keep HiGHS and preprocessing on the same default thread value.
        highs.setOptionValue("threads", preproc_thread_limit);
        highs_threads_option = static_cast<HighsInt>(preproc_thread_limit);
        if (threads_option_explicit) {
            // User provided a non-positive value (e.g. --threads 0): report the
            // resolved runtime value so settings output matches execution.
            accepted_highs_options["threads"] = std::to_string(preproc_thread_limit);
        } else {
            accepted_highs_options.erase("threads");
        }
    }
    preproc_thread_limit = std::max<int32_t>(1, preproc_thread_limit);
    if (threads_option_explicit) {
        accepted_highs_options["threads"] = std::to_string(preproc_thread_limit);
    }
    const unsigned logical_threads =
        std::max(1u, std::thread::hardware_concurrency());
    const unsigned physical_cores = detect_physical_cores();
    const char* tbb_str = ", TBB";
    const char* simd_str = "";
#ifdef __AVX512F__
    simd_str = ", AVX-512";
#elif defined(__AVX2__)
    simd_str = ", AVX2";
#elif defined(__AVX__)
    simd_str = ", AVX";
#elif defined(__SSE4_2__)
    simd_str = ", SSE4.2";
#endif
    logger_.log("Thread count: {} physical cores, {} logical processors, using up to {} threads{}{}",
                physical_cores, logical_threads, preproc_thread_limit, tbb_str, simd_str);

    std::vector<std::string> nondefault_settings;
    nondefault_settings.push_back("presolve=off");
    nondefault_settings.push_back("mip_pscost_minreliable=0");
    if (separation_interval != 1) {
        nondefault_settings.push_back(
            "separation_interval=" + std::to_string(separation_interval));
    }
    if (max_cuts_per_sep != 3) {
        nondefault_settings.push_back(
            "max_cuts_per_separator=" + std::to_string(max_cuts_per_sep));
    }
    if (all_pairs_propagation) nondefault_settings.push_back("all_pairs_propagation=on");
    if (disable_heuristics) nondefault_settings.push_back("disable_heuristics=on");
    if (!edge_elimination) nondefault_settings.push_back("edge_elimination=off");
    if (!edge_elimination_nodes) nondefault_settings.push_back("edge_elimination_nodes=off");
    if (branch_hyper != "off") nondefault_settings.push_back("branch_hyper=" + branch_hyper);
    if (deterministic_work_units > 0) {
        nondefault_settings.push_back(
            "deterministic_work_units=" + std::to_string(deterministic_work_units));
    }
    if (heuristic_node_interval != 200) {
        nondefault_settings.push_back(
            "heuristic_node_interval=" + std::to_string(heuristic_node_interval));
    }
    if (heuristic_deterministic_restarts != 32) {
        nondefault_settings.push_back(
            "heuristic_deterministic_restarts="
            + std::to_string(heuristic_deterministic_restarts));
    }
    if (preproc_stage1_bounds != "two_cycle") {
        nondefault_settings.push_back("preproc_stage1_bounds=" + preproc_stage1_bounds);
    }
    if (preproc_fast_restarts != 12) {
        nondefault_settings.push_back(
            "preproc_fast_restarts=" + std::to_string(preproc_fast_restarts));
    }
    {
        const preprocess::ng::DssrOptions defaults;
        if (dssr_options.initial_ng_size != defaults.initial_ng_size) {
            nondefault_settings.push_back(
                "ng_initial_size=" + std::to_string(dssr_options.initial_ng_size));
        }
        if (dssr_options.max_ng_size != defaults.max_ng_size) {
            nondefault_settings.push_back(
                "ng_max_size=" + std::to_string(dssr_options.max_ng_size));
        }
        if (dssr_options.dssr_iterations != defaults.dssr_iterations) {
            nondefault_settings.push_back(
                "ng_dssr_iters=" + std::to_string(dssr_options.dssr_iterations));
        }
        if (dssr_options.enable_simd != defaults.enable_simd) {
            nondefault_settings.push_back(
                std::string("ng_simd=") + (dssr_options.enable_simd ? "on" : "off"));
        }
    }
    if (cutoff < std::numeric_limits<double>::infinity()) {
        nondefault_settings.push_back("cutoff=" + std::to_string(cutoff));
    }
    for (const auto& [k, v] : accepted_highs_options) {
        nondefault_settings.push_back(k + "=" + v);
    }
    const std::string settings_summary = join_tokens(nondefault_settings);

    const int32_t source = problem_.source();
    const int32_t target = problem_.target();
    const int32_t n = problem_.num_nodes();
    const int32_t m = problem_.num_edges();
    const std::string solve_label =
        problem_.name.empty() ? std::string("instance") : problem_.name;
    logger_.log("Solving {} with:", solve_label);
    if (problem_.is_tour()) {
        logger_.log("  {} nodes, {} edges, tour through node {}", n, m, source);
    } else {
        logger_.log("  {} nodes, {} edges, path from node {} to node {}", n, m, source, target);
    }
    logger_.log("  Settings  : {}", settings_summary);

    auto non_highs_work_budget =
        std::make_shared<WorkUnitBudget>(deterministic_work_units);

    model_detail::SolveConcurrencyGuard solve_guard(max_concurrent_solves);

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

    // Honor user output preference for HiGHS logging (single-source output).
    highs.setOptionValue("output_flag", output_flag);

    // Deterministic staged preprocessing:
    // 1) Stage-1 bounds (parallel forward/backward)
    // 2) Construction seeds (1/2-customer)
    // 3) Local search on top seeds (parallel)
    // 3.5) Edge elimination from construction UB
    // 4) Reduced local search + deeper ng/DSSR (parallel)
    // 5) Start HiGHS
    // correction = profit(source) when s == t (tour: depot profit double-subtracted)
    double correction = problem_.is_tour() ? problem_.profit(source) : 0.0;

    std::vector<double> fwd_bounds, bwd_bounds;
    std::vector<double> all_pairs;
    heuristic::HeuristicResult warm_start{{}, std::numeric_limits<double>::infinity()};
    double warm_start_ub = std::numeric_limits<double>::infinity();
    int32_t ng_size_used = 1;
    double ng_path_ub = std::numeric_limits<double>::infinity();
    std::vector<int32_t> ng_path_nodes;
    int32_t stage3_elim_edges = 0;
    double stage3_elim_ratio = 0.0;
    int32_t stage5_elim_edges = 0;
    double stage5_elim_ratio = 0.0;

    // Use cutoff as UB when provided (known optimal for prove-only mode)
    if (cutoff < std::numeric_limits<double>::infinity()) {
        warm_start_ub = cutoff;
    }

    {
        Timer preproc_timer;
        enum class Stage1BoundsBackend {
            two_cycle,
            ng_dssr,
        };
        Stage1BoundsBackend stage1_backend = Stage1BoundsBackend::two_cycle;
        if (preproc_stage1_bounds == "ng_dssr") {
            stage1_backend = Stage1BoundsBackend::ng_dssr;
        } else if (preproc_stage1_bounds == "auto") {
            stage1_backend = ((dssr_options.max_ng_size > 1)
                              || (dssr_options.dssr_iterations > 1))
                ? Stage1BoundsBackend::ng_dssr
                : Stage1BoundsBackend::two_cycle;
        }
        const char* stage1_backend_name =
            (stage1_backend == Stage1BoundsBackend::ng_dssr) ? "ng_dssr" : "two_cycle";
        if (workflow_dump) {
            logger_.log("Workflow DAG (startup):");
            logger_.log("  Stage1 [Bounds] -> backend={}", stage1_backend_name);
            logger_.log("  Stage2 [Construction] -> seed 1/2-customer tours/paths");
            logger_.log("  Stage3 [EdgeElim on construction UB]");
            logger_.log("  Stage4 [Parallel local search on reduced graph] -> optional all-pairs bounds");
            logger_.log("  Stage5 [EdgeElim on local-search UB]");
            logger_.log("  Stage6 [HiGHS] -> highs.run()");
        }
        logger_.log("Startup Stage1 [Bounds]: backend={} (source={}, target={})",
                    stage1_backend_name, source, target);

        tbb::task_arena preproc_arena(preproc_thread_limit);
        preproc_arena.execute([&] {
            tbb::task_group stage1_tg;
            int32_t ng_size_source = 1;
            int32_t ng_size_target = 1;
            if (stage1_backend == Stage1BoundsBackend::two_cycle) {
                stage1_tg.run([&] {
                    fwd_bounds = preprocess::labeling_from(problem_, source);
                });
                if (source != target) {
                    stage1_tg.run([&] {
                        bwd_bounds = preprocess::labeling_from(problem_, target);
                    });
                }
            } else {
                auto stage1_opts = dssr_options;
                stage1_opts.initial_ng_size = std::max<int32_t>(0, dssr_options.initial_ng_size);
                stage1_opts.max_ng_size = std::max<int32_t>(
                    1, std::max<int32_t>(stage1_opts.initial_ng_size, dssr_options.max_ng_size));
                stage1_opts.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations);

                stage1_tg.run([&] {
                    auto source_bounds = preprocess::ng::compute_bounds(
                        problem_, source, source, stage1_opts);
                    fwd_bounds = std::move(source_bounds.fwd);
                    ng_size_source = source_bounds.ng_size;
                    if (source_bounds.elementary_path_found) {
                        ng_path_ub = source_bounds.elementary_path_cost;
                        ng_path_nodes = std::move(source_bounds.elementary_path);
                    }
                });
                if (source != target) {
                    stage1_tg.run([&] {
                        auto target_bounds = preprocess::ng::compute_bounds(
                            problem_, target, target, stage1_opts);
                        bwd_bounds = std::move(target_bounds.fwd);
                        ng_size_target = target_bounds.ng_size;
                    });
                }
            }
            stage1_tg.wait();
            ng_size_used = std::max<int32_t>(
                ng_size_used, std::max(ng_size_source, ng_size_target));
        });
        if (source == target) {
            bwd_bounds = fwd_bounds;
        }
        logger_.log(
            "Startup Stage1 [Bounds] done: fwd_finite={}/{}, bwd_finite={}/{}, fwd[target]={}",
            count_finite_bounds(fwd_bounds), n,
            count_finite_bounds(bwd_bounds), n,
            (target >= 0 && target < n && target < static_cast<int32_t>(fwd_bounds.size()))
                ? fwd_bounds[static_cast<size_t>(target)]
                : std::numeric_limits<double>::infinity());
        heuristic::ConstructionPool construction_pool;
        heuristic::HeuristicResult construction_start{{}, std::numeric_limits<double>::infinity()};
        heuristic::HeuristicResult ls_stage4{{}, std::numeric_limits<double>::infinity()};
        double construction_ub = std::numeric_limits<double>::infinity();
        std::vector<bool> stage3_edge_active;

        if (!disable_heuristics) {
            const int32_t candidate_cap = std::max<int32_t>(
                preproc_fast_restarts, preproc_thread_limit);
            construction_pool = heuristic::build_construction_pool(
                problem_, candidate_cap);
            construction_start = heuristic::best_construction_solution(problem_, construction_pool);
            if (has_feasible_heuristic(construction_start)) {
                construction_ub = construction_start.objective;
                warm_start = construction_start;
                logger_.log("Startup Stage2 [Construction] done: UB0={}, candidates={}",
                            construction_ub,
                            static_cast<int32_t>(construction_pool.candidates.size()));
            } else {
                logger_.log("Startup Stage2 [Construction] done: no feasible incumbent");
            }
        } else {
            logger_.log("Startup Stage2 [Construction] skipped: disable_heuristics=true");
        }

        if (cutoff >= std::numeric_limits<double>::infinity()) {
            warm_start_ub = has_feasible_heuristic(warm_start)
                ? warm_start.objective
                : std::numeric_limits<double>::infinity();
        } else if (has_feasible_heuristic(warm_start)) {
            warm_start_ub = std::min(warm_start_ub, warm_start.objective);
        }

        if (edge_elimination
            && std::isfinite(construction_ub)
            && !fwd_bounds.empty() && !bwd_bounds.empty()) {
            auto eliminated = preprocess::edge_elimination(
                problem_, fwd_bounds, bwd_bounds, construction_ub, correction);
            stage3_elim_edges = static_cast<int32_t>(
                std::count(eliminated.begin(), eliminated.end(), true));
            stage3_elim_ratio = (m > 0)
                ? static_cast<double>(stage3_elim_edges) / static_cast<double>(m)
                : 0.0;
            stage3_edge_active.assign(static_cast<size_t>(m), true);
            for (int32_t e = 0; e < m; ++e) {
                if (eliminated[static_cast<size_t>(e)]) {
                    stage3_edge_active[static_cast<size_t>(e)] = false;
                }
            }
            logger_.log(
                "Startup Stage3 [EdgeElim on construction UB]: eliminated={}/{} ({:.1f}%), ub0={}",
                stage3_elim_edges, m, 100.0 * stage3_elim_ratio, construction_ub);
        } else {
            logger_.log(
                "Startup Stage3 [EdgeElim on construction UB] skipped: edge_elimination={}, finite_ub0={}, have_bounds={}",
                edge_elimination ? "true" : "false",
                std::isfinite(construction_ub) ? "true" : "false",
                (!fwd_bounds.empty() && !bwd_bounds.empty()) ? "true" : "false");
        }

        preproc_arena.execute([&] {
            tbb::task_group stage4_tg;
            if (!disable_heuristics && !construction_pool.candidates.empty()) {
                logger_.log(
                    "Startup Stage4 [Parallel local search on reduced graph]: starts={}, max_iter=200",
                    preproc_thread_limit);
                stage4_tg.run([&] {
                    ls_stage4 = heuristic::run_local_search_from_pool(
                        problem_, construction_pool, preproc_thread_limit, 200,
                        non_highs_work_budget, preproc_thread_limit,
                        stage3_edge_active);
                });
            } else {
                logger_.log("Startup Stage4 [Parallel local search on reduced graph] skipped");
            }

            if (all_pairs_propagation) {
                logger_.log("Startup Stage4 [Parallel local search on reduced graph]: all-pairs 2-cycle bounds running");
                stage4_tg.run([&] {
                    constexpr double inf = std::numeric_limits<double>::infinity();
                    all_pairs.assign(static_cast<size_t>(n) * n, inf);
                    tbb::parallel_for(0, n, [&](int32_t s) {
                        auto row = preprocess::forward_labeling(problem_, s);
                        std::copy(row.begin(), row.end(),
                                  all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
                    });
                    logger_.log("Startup Stage4 [Parallel local search on reduced graph]: all-pairs 2-cycle bounds done");
                });
            }
            stage4_tg.wait();
        });

        if (has_feasible_heuristic(ls_stage4) &&
            (!has_feasible_heuristic(warm_start) ||
             ls_stage4.objective + 1e-9 < warm_start.objective)) {
            logger_.log("Startup Stage4 [Parallel local search on reduced graph] improved UB: {} -> {}",
                        warm_start.objective, ls_stage4.objective);
            warm_start = std::move(ls_stage4);
        } else if (has_feasible_heuristic(ls_stage4)) {
            logger_.log("Startup Stage4 [Parallel local search on reduced graph] done: objective={} (no UB improvement)",
                        ls_stage4.objective);
        }

        if (all_pairs_propagation && !all_pairs.empty()) {
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
        }

        if (cutoff >= std::numeric_limits<double>::infinity()) {
            warm_start_ub = has_feasible_heuristic(warm_start)
                ? warm_start.objective
                : std::numeric_limits<double>::infinity();
        } else if (has_feasible_heuristic(warm_start)) {
            warm_start_ub = std::min(warm_start_ub, warm_start.objective);
        }
        if (ng_path_ub < warm_start_ub) {
            warm_start_ub = ng_path_ub;
        }
        if (!ng_path_nodes.empty()) {
            auto ng_start = build_ng_path_start(problem_, ng_path_nodes);
            if (has_feasible_heuristic(ng_start) &&
                (!has_feasible_heuristic(warm_start) ||
                 ng_start.objective < warm_start.objective - 1e-9)) {
                warm_start = std::move(ng_start);
                warm_start_ub = std::min(warm_start_ub, warm_start.objective);
            }
        }

        if (edge_elimination
            && std::isfinite(warm_start_ub)
            && !fwd_bounds.empty() && !bwd_bounds.empty()) {
            auto eliminated = preprocess::edge_elimination(
                problem_, fwd_bounds, bwd_bounds, warm_start_ub, correction);
            stage5_elim_edges = static_cast<int32_t>(
                std::count(eliminated.begin(), eliminated.end(), true));
            stage5_elim_ratio = (m > 0)
                ? static_cast<double>(stage5_elim_edges) / static_cast<double>(m)
                : 0.0;
            logger_.log(
                "Startup Stage5 [EdgeElim on local-search UB]: eliminated={}/{} ({:.1f}%), ub={}",
                stage5_elim_edges, m, 100.0 * stage5_elim_ratio, warm_start_ub);
        } else {
            logger_.log(
                "Startup Stage5 [EdgeElim on local-search UB] skipped: edge_elimination={}, finite_ub={}, have_bounds={}",
                edge_elimination ? "true" : "false",
                std::isfinite(warm_start_ub) ? "true" : "false",
                (!fwd_bounds.empty() && !bwd_bounds.empty()) ? "true" : "false");
        }

        logger_.log("Preprocessing: {}s, UB={}{} [Stage1 [Bounds]={}, Stage3 [EdgeElim on construction UB]: {}/{} ({:.1f}%), Stage5 [EdgeElim on local-search UB]: {}/{} ({:.1f}%)]",
                    preproc_timer.elapsed_seconds(),
                    warm_start_ub,
                    (cutoff < std::numeric_limits<double>::infinity() ? " (cutoff)" : ""),
                    stage1_backend_name,
                    stage3_elim_edges, m, 100.0 * stage3_elim_ratio,
                    stage5_elim_edges, m, 100.0 * stage5_elim_ratio);
    }
    std::vector<int8_t> fixed_y(static_cast<size_t>(n), static_cast<int8_t>(-1));
    HiGHSBridge bridge(problem_, highs, logger_, separation_tol);
    bridge.set_separation_interval(separation_interval);
    bridge.set_max_cuts_per_separator(max_cuts_per_sep);
    bridge.set_submip_separation(submip_separation);
    bridge.set_upper_bound(warm_start_ub);
    bridge.set_rc_fixing(rc_settings);
    bridge.set_edge_elimination(edge_elimination);
    bridge.set_edge_elimination_nodes(edge_elimination_nodes);
    bridge.set_heuristic_node_interval(heuristic_node_interval);
    bridge.set_heuristic_deterministic_restarts(heuristic_deterministic_restarts);
    bridge.set_fixed_y(fixed_y);
    bridge.set_work_unit_budget(non_highs_work_budget);
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
    bridge.set_shared_bounds_store(shared_bounds);
    if (workflow_dump) {
        logger_.log("Workflow DAG (solve):");
        logger_.log("  build_formulation -> highs.run() with callbacks");
    }

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
    {
        const auto& lp = highs.getLp();
        constexpr double tol = 1e-9;
        int64_t n_binary = 0;
        int64_t n_integer = 0;
        int64_t n_continuous = 0;
        for (HighsInt j = 0; j < lp.num_col_; ++j) {
            HighsVarType vtype = HighsVarType::kContinuous;
            if (j >= 0 && j < static_cast<HighsInt>(lp.integrality_.size())) {
                vtype = lp.integrality_[j];
            }
            const bool integral = (vtype == HighsVarType::kInteger
                                   || vtype == HighsVarType::kSemiInteger);
            if (!integral) {
                n_continuous++;
                continue;
            }
            const double lb = lp.col_lower_[j];
            const double ub = lp.col_upper_[j];
            const bool is_binary = lb >= -tol && lb <= 1.0 + tol
                && ub >= -tol && ub <= 1.0 + tol;
            if (is_binary) n_binary++;
            else n_integer++;
        }
        std::vector<std::string> type_tokens;
        if (n_binary > 0) type_tokens.push_back(std::to_string(n_binary) + " binary");
        if (n_integer > 0) type_tokens.push_back(std::to_string(n_integer) + " integer");
        if (n_continuous > 0) type_tokens.push_back(std::to_string(n_continuous) + " continuous");
        const std::string type_summary =
            type_tokens.empty() ? "" : (" (" + join_tokens(type_tokens, ", ") + ")");
        logger_.log("Built MIP model:");
        logger_.log("  {} rows, {} cols{}, {} nonzeros",
                    lp.num_row_, lp.num_col_, type_summary, lp.a_matrix_.numNz());
    }

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
    } else {
        // Ensure no stale global hyperplane callback leaks from prior solves.
        HighsUserSeparator::setBranchingCallback({});
        HighsUserSeparator::clearStack();
    }

    bridge.install_separators();
    bridge.install_propagator();
    bridge.set_heuristic_callback(heuristic_callback);
    bridge.set_heuristic_budget_ms(heuristic_budget_ms);
    bridge.set_heuristic_strategy(heuristic_strategy);
    bridge.install_heuristic_callback();

    if (!disable_heuristics && has_feasible_heuristic(warm_start)) {
        // Pass heuristic solution to HiGHS (skip if heuristics disabled)
        HighsSolution start;
        start.value_valid = true;
        start.col_value = std::move(warm_start.col_values);
        HighsInt num_cols = highs.getNumCol();
        while (static_cast<HighsInt>(start.col_value.size()) < num_cols)
            start.col_value.push_back(0.0);
        highs.setSolution(start);
    }

    logger_.log("Startup Stage6 [HiGHS]: launching highs.run()");
    highs.run();

    // Clear static callbacks/state to avoid leaking between solves
    HighsUserSeparator::clearCallback();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();
    if (non_highs_work_budget && non_highs_work_budget->used() > 0) {
        if (non_highs_work_budget->capped()) {
            logger_.log("Non-HiGHS work units used: {}/{}",
                        non_highs_work_budget->used(),
                        non_highs_work_budget->limit());
        } else {
            logger_.log("Non-HiGHS work units used: {} (uncapped)",
                        non_highs_work_budget->used());
        }
    }

    return result;
}

}  // namespace rcspp
