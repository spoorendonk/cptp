#include "model/model.h"
#include "model/deterministic_checkpoint.h"
#include "model/highs_bridge.h"
#include "model/solve_concurrency_guard.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <thread>
#include <unordered_map>

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

constexpr uint64_t kDssrSigOffset = 1469598103934665603ull;
constexpr uint64_t kDssrSigPrime = 1099511628211ull;

uint64_t dssr_sig_mix(uint64_t sig, uint64_t value) {
    sig ^= value;
    sig *= kDssrSigPrime;
    return sig;
}

uint64_t dssr_sig_quantize(double value) {
    if (!std::isfinite(value)) {
        return std::numeric_limits<uint64_t>::max();
    }
    return std::bit_cast<uint64_t>(value);
}

struct ParamipFixing {
    int32_t node = -1;
    bool value = false;
};

struct ParamipChunkPlan {
    int32_t id = 0;
    std::vector<ParamipFixing> fixings;
};

int32_t ceil_log2_chunks(int32_t chunks) {
    int32_t depth = 0;
    int32_t cap = 1;
    while (cap < chunks && depth < 30) {
        cap <<= 1;
        ++depth;
    }
    return depth;
}

std::vector<int32_t> choose_paramip_split_nodes(const Problem& problem,
                                                int32_t depth,
                                                std::span<const int8_t> fixed_y = {}) {
    std::vector<int32_t> candidates;
    const int32_t n = problem.num_nodes();
    const int32_t source = problem.source();
    const int32_t target = problem.target();
    candidates.reserve(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) {
        if (i == source || i == target) continue;
        if (!fixed_y.empty() &&
            i < static_cast<int32_t>(fixed_y.size()) &&
            fixed_y[static_cast<size_t>(i)] != static_cast<int8_t>(-1)) {
            continue;
        }
        candidates.push_back(i);
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [&](int32_t a, int32_t b) {
                         const double pa = problem.profit(a);
                         const double pb = problem.profit(b);
                         if (pa != pb) return pa > pb;
                         return a < b;
                     });
    if (depth < static_cast<int32_t>(candidates.size())) {
        candidates.resize(static_cast<size_t>(depth));
    }
    return candidates;
}

std::vector<ParamipChunkPlan> build_paramip_chunk_plan(
    std::span<const int32_t> split_nodes) {
    const int32_t depth = static_cast<int32_t>(split_nodes.size());
    std::vector<ParamipChunkPlan> chunks;
    if (depth <= 0) return chunks;
    const int32_t total = 1 << depth;
    chunks.reserve(static_cast<size_t>(total));
    for (int32_t cid = 0; cid < total; ++cid) {
        ParamipChunkPlan chunk;
        chunk.id = cid;
        chunk.fixings.reserve(static_cast<size_t>(depth));
        for (int32_t level = 0; level < depth; ++level) {
            const bool value = ((cid >> level) & 1) != 0;
            chunk.fixings.push_back(ParamipFixing{
                .node = split_nodes[static_cast<size_t>(level)],
                .value = value,
            });
        }
        chunks.push_back(std::move(chunk));
    }
    return chunks;
}

std::string encode_paramip_fixings(std::span<const int8_t> fixed_y) {
    std::ostringstream ss;
    bool first = true;
    for (size_t i = 0; i < fixed_y.size(); ++i) {
        const int8_t v = fixed_y[i];
        if (v != static_cast<int8_t>(0) && v != static_cast<int8_t>(1)) continue;
        if (!first) ss << ",";
        first = false;
        ss << static_cast<int32_t>(i) << ":" << static_cast<int32_t>(v);
    }
    return ss.str();
}

std::vector<int8_t> parse_paramip_fixings(const std::string& spec,
                                          int32_t num_nodes,
                                          int32_t source,
                                          int32_t target,
                                          Logger& logger) {
    std::vector<int8_t> fixed_y(static_cast<size_t>(num_nodes), static_cast<int8_t>(-1));
    if (spec.empty()) return fixed_y;
    std::stringstream ss(spec);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        const auto pos = token.find(':');
        if (pos == std::string::npos) {
            logger.log("Warning: invalid paramip_fixings token '{}', expected node:value", token);
            continue;
        }
        const std::string node_str = token.substr(0, pos);
        const std::string value_str = token.substr(pos + 1);
        int node = -1;
        int value = -1;
        try {
            node = std::stoi(node_str);
            value = std::stoi(value_str);
        } catch (...) {
            logger.log("Warning: invalid paramip_fixings token '{}', expected integers", token);
            continue;
        }
        if (node < 0 || node >= num_nodes || (value != 0 && value != 1)) {
            logger.log("Warning: invalid paramip_fixings token '{}', out of range", token);
            continue;
        }
        if (node == source || node == target) {
            if (value == 0) {
                logger.log("Warning: ignoring paramip_fixings on terminal node {}:{}", node, value);
                continue;
            }
            fixed_y[static_cast<size_t>(node)] = 1;
            continue;
        }
        fixed_y[static_cast<size_t>(node)] = static_cast<int8_t>(value);
    }
    return fixed_y;
}

SolverOptions make_paramip_subsolve_options(const SolverOptions& base_options,
                                            std::string fixings_spec,
                                            double cutoff,
                                            double remaining_time_limit_sec,
                                            bool force_quiet,
                                            std::optional<int32_t> random_seed = std::nullopt,
                                            std::optional<int32_t> mip_max_nodes = std::nullopt) {
    SolverOptions out;
    out.reserve(base_options.size() + 12);
    for (const auto& [k, v] : base_options) {
        if (k == "paramip_mode" || k == "paramip_chunks" || k == "paramip_workers" ||
            k == "paramip_fixings" || k == "cutoff" || k == "threads" ||
            k == "time_limit" ||
            k == "max_concurrent_solves" || k == "branch_hyper" ||
            k == "paramip_root_probes" || k == "paramip_root_pick" ||
            k == "_internal_skip_preproc" ||
            k == "_internal_prepared_ctx_id" ||
            k == "_internal_disable_async_dssr" ||
            (k == "random_seed" && random_seed.has_value()) ||
            (k == "mip_max_nodes" && mip_max_nodes.has_value())) {
            continue;
        }
        if (force_quiet && k == "output_flag") continue;
        out.push_back({k, v});
    }
    out.push_back({"paramip_mode", "off"});
    out.push_back({"paramip_fixings", std::move(fixings_spec)});
    out.push_back({"threads", "1"});
    out.push_back({"max_concurrent_solves", "0"});
    out.push_back({"branch_hyper", "off"});
    if (force_quiet) {
        out.push_back({"output_flag", "false"});
    }
    if (std::isfinite(cutoff)) {
        out.push_back({"cutoff", std::to_string(cutoff)});
    }
    if (std::isfinite(remaining_time_limit_sec)) {
        out.push_back({"time_limit", std::to_string(std::max(0.001, remaining_time_limit_sec))});
    }
    if (random_seed.has_value()) {
        out.push_back({"random_seed", std::to_string(*random_seed)});
    }
    if (mip_max_nodes.has_value()) {
        out.push_back({"mip_max_nodes", std::to_string(*mip_max_nodes)});
    }
    return out;
}

void aggregate_separator_stats(std::map<std::string, SeparatorStats>& dst,
                               const std::map<std::string, SeparatorStats>& src) {
    for (const auto& [name, s] : src) {
        auto& d = dst[name];
        d.cuts_added += s.cuts_added;
        d.rounds_called += s.rounds_called;
        d.time_seconds += s.time_seconds;
    }
}

struct PreparedSolveContext {
    std::vector<double> fwd_bounds;
    std::vector<double> bwd_bounds;
    std::vector<double> all_pairs;
    heuristic::HeuristicResult warm_start;
    double warm_start_ub = std::numeric_limits<double>::infinity();
    int32_t ng_size_used = 1;
    double correction = 0.0;
    int32_t stage1_elim_edges = 0;
    double stage1_elim_ratio = 0.0;
    std::string stage1_backend_name = "reused";
};

std::mutex g_prepared_context_mutex;
std::unordered_map<int64_t, std::shared_ptr<const PreparedSolveContext>>
    g_prepared_contexts;
std::atomic<int64_t> g_next_prepared_context_id{1};

int64_t register_prepared_context(std::shared_ptr<const PreparedSolveContext> ctx) {
    const int64_t id = g_next_prepared_context_id.fetch_add(1, std::memory_order_relaxed);
    std::lock_guard<std::mutex> lock(g_prepared_context_mutex);
    g_prepared_contexts[id] = std::move(ctx);
    return id;
}

std::shared_ptr<const PreparedSolveContext> find_prepared_context(int64_t id) {
    if (id <= 0) return {};
    std::lock_guard<std::mutex> lock(g_prepared_context_mutex);
    auto it = g_prepared_contexts.find(id);
    if (it == g_prepared_contexts.end()) return {};
    return it->second;
}

void unregister_prepared_context(int64_t id) {
    if (id <= 0) return;
    std::lock_guard<std::mutex> lock(g_prepared_context_mutex);
    g_prepared_contexts.erase(id);
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
    std::string parallel_mode = "deterministic";
    bool deterministic = true;
    int64_t deterministic_work_units = 0;
    int32_t heuristic_deterministic_restarts = 32;
    int64_t heuristic_node_interval = 200;
    bool heuristic_async_injection = true;
    bool workflow_dump = false;
    std::string paramip_mode = "off";
    int32_t paramip_chunks = 0;
    int32_t paramip_workers = 0;
    int32_t paramip_root_probes = 0;
    std::string paramip_root_pick = "auto";
    std::string paramip_fixings_spec;
    bool internal_skip_preproc = false;
    int64_t internal_prepared_ctx_id = 0;
    bool internal_disable_async_dssr = false;
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
    bool dssr_background_updates = true;
    bool dssr_background_updates_explicit = false;
    std::string dssr_background_policy = "fixed";
    int32_t dssr_background_max_epochs = 0;
    int32_t dssr_background_auto_min_epochs = 4;
    int32_t dssr_background_auto_no_progress_limit = 6;
    bool preproc_adaptive = true;
    std::string preproc_stage1_bounds = "auto";
    int32_t preproc_fast_restarts = 12;
    double preproc_fast_budget_ms = 30.0;
    int32_t preproc_second_ws_large_n = 60;
    double preproc_second_ws_min_elim = 0.05;
    double preproc_second_ws_min_elim_large = 0.02;
    double preproc_second_ws_budget_ms_min = 20.0;
    double preproc_second_ws_budget_ms_max = 400.0;
    double preproc_second_ws_budget_scale = 8.0;
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
        if (key == "workflow_dump") {
            workflow_dump = (value == "true" || value == "1");
            continue;
        }
        if (key == "paramip_mode") {
            std::string mode = value;
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (mode == "off" || mode == "plan" || mode == "static_root") {
                paramip_mode = mode;
            } else {
                logger_.log("Warning: unknown paramip_mode='{}' (expected off|plan|static_root), using off",
                            value);
                paramip_mode = "off";
            }
            continue;
        }
        if (key == "paramip_chunks") {
            paramip_chunks = std::stoi(value);
            continue;
        }
        if (key == "paramip_workers") {
            paramip_workers = std::stoi(value);
            continue;
        }
        if (key == "paramip_root_probes") {
            paramip_root_probes = std::stoi(value);
            continue;
        }
        if (key == "paramip_root_pick") {
            std::string pick = value;
            std::transform(pick.begin(), pick.end(), pick.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (pick == "auto" || pick == "best" || pick == "first") {
                paramip_root_pick = pick;
            } else {
                logger_.log("Warning: unknown paramip_root_pick='{}' (expected auto|best|first), using auto",
                            value);
                paramip_root_pick = "auto";
            }
            continue;
        }
        if (key == "paramip_fixings") {
            paramip_fixings_spec = value;
            continue;
        }
        if (key == "_internal_skip_preproc") {
            internal_skip_preproc = (value == "true" || value == "1");
            continue;
        }
        if (key == "_internal_prepared_ctx_id") {
            try {
                internal_prepared_ctx_id = std::stoll(value);
            } catch (...) {
                internal_prepared_ctx_id = 0;
            }
            continue;
        }
        if (key == "_internal_disable_async_dssr") {
            internal_disable_async_dssr = (value == "true" || value == "1");
            continue;
        }
        if (key == "parallel_mode") {
            std::string mode = value;
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (mode == "deterministic" || mode == "opportunistic") {
                parallel_mode = mode;
            } else {
                logger_.log("Warning: unknown parallel_mode='{}' (expected deterministic|opportunistic), using deterministic",
                            value);
                parallel_mode = "deterministic";
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
        if (key == "dssr_background_updates") {
            dssr_background_updates = (value == "true" || value == "1");
            dssr_background_updates_explicit = true;
            continue;
        }
        if (key == "dssr_background_policy") {
            std::string policy = value;
            std::transform(policy.begin(), policy.end(), policy.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (policy == "fixed" || policy == "auto") {
                dssr_background_policy = policy;
            } else {
                logger_.log("Warning: unknown dssr_background_policy='{}' (expected fixed|auto), using fixed",
                            value);
                dssr_background_policy = "fixed";
            }
            continue;
        }
        if (key == "dssr_background_max_epochs") {
            dssr_background_max_epochs = std::stoi(value);
            continue;
        }
        if (key == "dssr_background_auto_min_epochs") {
            dssr_background_auto_min_epochs = std::stoi(value);
            continue;
        }
        if (key == "dssr_background_auto_no_progress_limit") {
            dssr_background_auto_no_progress_limit = std::stoi(value);
            continue;
        }
        if (key == "preproc_adaptive") {
            preproc_adaptive = (value == "true" || value == "1");
            continue;
        }
        if (key == "preproc_stage1_bounds") {
            std::string mode = value;
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (mode == "two_cycle" || mode == "ng1" || mode == "auto") {
                preproc_stage1_bounds = mode;
            } else {
                logger_.log("Warning: unknown preproc_stage1_bounds='{}' (expected two_cycle|ng1|auto), using auto",
                            value);
                preproc_stage1_bounds = "auto";
            }
            continue;
        }
        if (key == "preproc_fast_restarts") {
            preproc_fast_restarts = std::stoi(value);
            continue;
        }
        if (key == "preproc_fast_budget_ms") {
            preproc_fast_budget_ms = std::stod(value);
            continue;
        }
        if (key == "preproc_second_ws_large_n") {
            preproc_second_ws_large_n = std::stoi(value);
            continue;
        }
        if (key == "preproc_second_ws_min_elim") {
            preproc_second_ws_min_elim = std::stod(value);
            continue;
        }
        if (key == "preproc_second_ws_min_elim_large") {
            preproc_second_ws_min_elim_large = std::stod(value);
            continue;
        }
        if (key == "preproc_second_ws_budget_ms_min") {
            preproc_second_ws_budget_ms_min = std::stod(value);
            continue;
        }
        if (key == "preproc_second_ws_budget_ms_max") {
            preproc_second_ws_budget_ms_max = std::stod(value);
            continue;
        }
        if (key == "preproc_second_ws_budget_scale") {
            preproc_second_ws_budget_scale = std::stod(value);
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
        }
    }

    dssr_options.initial_ng_size = std::max<int32_t>(1, dssr_options.initial_ng_size);
    dssr_options.max_ng_size = std::max<int32_t>(dssr_options.initial_ng_size,
                                                 dssr_options.max_ng_size);
    dssr_options.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations);
    preproc_fast_restarts = std::max<int32_t>(1, preproc_fast_restarts);
    preproc_fast_budget_ms = std::max(0.0, preproc_fast_budget_ms);
    preproc_second_ws_large_n = std::max<int32_t>(1, preproc_second_ws_large_n);
    preproc_second_ws_min_elim = std::clamp(preproc_second_ws_min_elim, 0.0, 1.0);
    preproc_second_ws_min_elim_large = std::clamp(preproc_second_ws_min_elim_large, 0.0, 1.0);
    preproc_second_ws_budget_ms_min = std::max(0.0, preproc_second_ws_budget_ms_min);
    preproc_second_ws_budget_ms_max = std::max(preproc_second_ws_budget_ms_min,
                                                preproc_second_ws_budget_ms_max);
    preproc_second_ws_budget_scale = std::max(0.0, preproc_second_ws_budget_scale);
    heuristic_node_interval = std::max<int64_t>(1, heuristic_node_interval);
    deterministic = (parallel_mode != "opportunistic");
    deterministic_work_units = std::max<int64_t>(0, deterministic_work_units);
    heuristic_deterministic_restarts =
        std::max<int32_t>(1, heuristic_deterministic_restarts);
    max_concurrent_solves = std::max<int32_t>(0, max_concurrent_solves);
    paramip_chunks = std::max<int32_t>(0, paramip_chunks);
    paramip_workers = std::max<int32_t>(0, paramip_workers);
    paramip_root_probes = std::max<int32_t>(0, paramip_root_probes);
    hyper_sb_iter_limit = std::max<int32_t>(1, hyper_sb_iter_limit);
    hyper_sb_min_reliable = std::max<int32_t>(1, hyper_sb_min_reliable);
    hyper_sb_max_candidates = std::max<int32_t>(1, hyper_sb_max_candidates);
    hyper_sb_max_depth = std::max<int32_t>(0, hyper_sb_max_depth);
    dssr_background_max_epochs = std::max<int32_t>(0, dssr_background_max_epochs);
    dssr_background_auto_min_epochs =
        std::max<int32_t>(1, dssr_background_auto_min_epochs);
    dssr_background_auto_no_progress_limit =
        std::max<int32_t>(1, dssr_background_auto_no_progress_limit);
    if (deterministic && !dssr_background_updates_explicit) {
        dssr_background_updates = false;
    }

    const int32_t source = problem_.source();
    const int32_t target = problem_.target();
    const int32_t n = problem_.num_nodes();

    auto fixed_y = parse_paramip_fixings(
        paramip_fixings_spec, n, source, target, logger_);
    std::shared_ptr<const PreparedSolveContext> prepared_context;
    if (internal_skip_preproc) {
        prepared_context = find_prepared_context(internal_prepared_ctx_id);
        if (!prepared_context) {
            logger_.log(
                "Warning: internal prepared context {} not found; running normal preprocessing",
                internal_prepared_ctx_id);
            internal_skip_preproc = false;
        }
    }

    int32_t base_random_seed = 0;
    for (const auto& [k, v] : options) {
        if (k != "random_seed") continue;
        try {
            base_random_seed = std::stoi(v);
        } catch (...) {
            base_random_seed = 0;
        }
        break;
    }

    const int32_t paramip_worker_target =
        (paramip_workers > 0)
            ? paramip_workers
            : std::max<int32_t>(1, static_cast<int32_t>(std::thread::hardware_concurrency()));

    std::vector<ParamipChunkPlan> chunk_plan;
    if (paramip_mode != "off") {
        const int32_t requested_chunks = std::max<int32_t>(2, paramip_chunks);
        const int32_t depth = ceil_log2_chunks(requested_chunks);
        auto split_nodes = choose_paramip_split_nodes(problem_, depth, fixed_y);
        chunk_plan = build_paramip_chunk_plan(split_nodes);
        logger_.log("ParaMIP plan: mode={}, requested_chunks={}, depth={}, chunks={}, workers={}",
                    paramip_mode, requested_chunks, depth,
                    static_cast<int32_t>(chunk_plan.size()), paramip_worker_target);
        if (paramip_mode == "static_root") {
            logger_.log("ParaMIP stage0 config: root_probes={}, pick={}",
                        paramip_root_probes, paramip_root_pick);
        }
        if (static_cast<int32_t>(split_nodes.size()) < depth) {
            logger_.log("ParaMIP plan: limited split nodes ({}/{}) due to instance size",
                        static_cast<int32_t>(split_nodes.size()), depth);
        }
        if (workflow_dump && !chunk_plan.empty()) {
            const int32_t dump_limit = std::min<int32_t>(
                8, static_cast<int32_t>(chunk_plan.size()));
            for (int32_t i = 0; i < dump_limit; ++i) {
                const auto& chunk = chunk_plan[static_cast<size_t>(i)];
                std::ostringstream ss;
                ss << "  chunk#" << chunk.id << ": ";
                for (size_t k = 0; k < chunk.fixings.size(); ++k) {
                    if (k > 0) ss << ", ";
                    ss << "y[" << chunk.fixings[k].node << "]="
                       << (chunk.fixings[k].value ? 1 : 0);
                }
                logger_.log("{}", ss.str());
            }
        }
    }

    if (paramip_mode == "static_root" && !chunk_plan.empty() && internal_skip_preproc) {
        Timer paramip_timer;
        std::vector<SolveResult> chunk_results(chunk_plan.size());
        std::atomic<double> global_ub{cutoff};
        double global_time_limit_sec = std::numeric_limits<double>::infinity();
        for (const auto& [k, v] : options) {
            if (k != "time_limit") continue;
            try {
                global_time_limit_sec = std::stod(v);
            } catch (...) {
                global_time_limit_sec = std::numeric_limits<double>::infinity();
            }
            break;
        }

        int32_t selected_root_seed = base_random_seed;
        const std::string root_pick_policy =
            (paramip_root_pick == "auto")
                ? (deterministic ? "best" : "first")
                : paramip_root_pick;
        const int32_t root_probe_count =
            (paramip_root_probes > 0)
                ? paramip_root_probes
                : std::min<int32_t>(std::max<int32_t>(2, paramip_worker_target), 8);

        if (root_probe_count > 1) {
            struct RootProbeResult {
                int32_t probe_idx = -1;
                int32_t seed = 0;
                int32_t completion_rank = std::numeric_limits<int32_t>::max();
                SolveResult result;
            };

            std::vector<RootProbeResult> probes(static_cast<size_t>(root_probe_count));
            for (int32_t i = 0; i < root_probe_count; ++i) {
                probes[static_cast<size_t>(i)].probe_idx = i;
                probes[static_cast<size_t>(i)].seed = base_random_seed + i;
            }

            double stage0_remaining = std::numeric_limits<double>::infinity();
            if (std::isfinite(global_time_limit_sec)) {
                stage0_remaining = global_time_limit_sec - paramip_timer.elapsed_seconds();
            }
            if (stage0_remaining > 1e-6 || !std::isfinite(stage0_remaining)) {
                const double stage0_probe_limit =
                    std::isfinite(stage0_remaining)
                        ? std::max(0.05, std::min(5.0, 0.25 * stage0_remaining))
                        : std::numeric_limits<double>::infinity();
                const std::string root_fixings = encode_paramip_fixings(fixed_y);
                std::atomic<int32_t> completion_counter{0};
                auto run_probe = [&](int32_t idx) {
                    auto& probe = probes[static_cast<size_t>(idx)];
                    auto sub_opts = make_paramip_subsolve_options(
                        options, root_fixings, global_ub.load(std::memory_order_relaxed),
                        stage0_probe_limit, true,
                        probe.seed, 1);
                    Model sub_model;
                    sub_model.set_problem(problem_);
                    probe.result = sub_model.solve(sub_opts);
                    probe.completion_rank = completion_counter.fetch_add(
                        1, std::memory_order_relaxed);
                };

                const int32_t workers = std::min<int32_t>(
                    paramip_worker_target, root_probe_count);
                if (workers <= 1) {
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        run_probe(i);
                    }
                } else {
                    std::atomic<int32_t> next{0};
                    std::vector<std::thread> pool;
                    pool.reserve(static_cast<size_t>(workers));
                    for (int32_t w = 0; w < workers; ++w) {
                        pool.emplace_back([&] {
                            for (;;) {
                                const int32_t idx =
                                    next.fetch_add(1, std::memory_order_relaxed);
                                if (idx >= root_probe_count) break;
                                run_probe(idx);
                            }
                        });
                    }
                    for (auto& t : pool) t.join();
                }

                auto has_finite_bound = [](const RootProbeResult& p) {
                    return std::isfinite(p.result.bound);
                };
                auto has_incumbent = [](const RootProbeResult& p) {
                    if (!std::isfinite(p.result.objective)) return false;
                    if (p.result.status == SolveResult::Status::Optimal ||
                        p.result.status == SolveResult::Status::Feasible) {
                        return true;
                    }
                    if (p.result.status == SolveResult::Status::TimeLimit) {
                        return !p.result.tour.empty() || !p.result.tour_arcs.empty();
                    }
                    return false;
                };

                int32_t chosen = -1;
                if (root_pick_policy == "first") {
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        const auto& p = probes[static_cast<size_t>(i)];
                        if (p.result.status == SolveResult::Status::Error) continue;
                        if (!has_finite_bound(p)) continue;
                        if (chosen < 0 ||
                            p.completion_rank <
                                probes[static_cast<size_t>(chosen)].completion_rank) {
                            chosen = i;
                        }
                    }
                    if (chosen < 0) {
                        for (int32_t i = 0; i < root_probe_count; ++i) {
                            const auto& p = probes[static_cast<size_t>(i)];
                            if (p.result.status == SolveResult::Status::Error) continue;
                            if (chosen < 0 ||
                                p.completion_rank <
                                    probes[static_cast<size_t>(chosen)].completion_rank) {
                                chosen = i;
                            }
                        }
                    }
                } else {
                    // "best": pick strongest root LP bound (max bound in minimization),
                    // tie-breaking by incumbent objective then probe index.
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        const auto& p = probes[static_cast<size_t>(i)];
                        if (p.result.status == SolveResult::Status::Error) continue;
                        if (chosen < 0) {
                            chosen = i;
                            continue;
                        }
                        const auto& best = probes[static_cast<size_t>(chosen)];
                        const bool p_has_bound = has_finite_bound(p);
                        const bool b_has_bound = has_finite_bound(best);
                        if (p_has_bound && !b_has_bound) {
                            chosen = i;
                            continue;
                        }
                        if (p_has_bound && b_has_bound &&
                            p.result.bound > best.result.bound + 1e-9) {
                            chosen = i;
                            continue;
                        }
                        if (p_has_bound && b_has_bound &&
                            std::abs(p.result.bound - best.result.bound) <= 1e-9) {
                            const bool p_has_sol = has_incumbent(p);
                            const bool b_has_sol = has_incumbent(best);
                            if (p_has_sol && !b_has_sol) {
                                chosen = i;
                                continue;
                            }
                            if (p_has_sol && b_has_sol &&
                                p.result.objective < best.result.objective - 1e-9) {
                                chosen = i;
                                continue;
                            }
                        }
                    }
                }

                if (chosen >= 0) {
                    const auto& picked = probes[static_cast<size_t>(chosen)];
                    selected_root_seed = picked.seed;
                    double best_probe_obj = std::numeric_limits<double>::infinity();
                    for (const auto& p : probes) {
                        if (!has_incumbent(p)) continue;
                        best_probe_obj = std::min(best_probe_obj, p.result.objective);
                    }
                    if (std::isfinite(best_probe_obj)) {
                        double prev = global_ub.load(std::memory_order_relaxed);
                        while (best_probe_obj + 1e-9 < prev) {
                            if (global_ub.compare_exchange_weak(
                                    prev, best_probe_obj,
                                    std::memory_order_relaxed,
                                    std::memory_order_relaxed)) {
                                break;
                            }
                        }
                    }
                    logger_.log(
                        "ParaMIP stage0: probes={}, policy={}, chosen_probe={}, seed={}, bound={}, obj={}, status={}",
                        root_probe_count, root_pick_policy, chosen, selected_root_seed,
                        picked.result.bound, picked.result.objective,
                        static_cast<int>(picked.result.status));
                } else {
                    logger_.log(
                        "ParaMIP stage0: probes={}, policy={}, no valid root candidate (keeping seed={})",
                        root_probe_count, root_pick_policy, selected_root_seed);
                }
            }
        }

        auto encode_chunk_fixings = [&](const ParamipChunkPlan& chunk) {
            std::vector<int8_t> merged = fixed_y;
            if (merged.empty()) {
                merged.assign(static_cast<size_t>(n), static_cast<int8_t>(-1));
            }
            for (const auto& f : chunk.fixings) {
                if (f.node < 0 || f.node >= n) continue;
                merged[static_cast<size_t>(f.node)] = static_cast<int8_t>(f.value ? 1 : 0);
            }
            return encode_paramip_fixings(merged);
        };

        auto run_chunk = [&](int32_t idx) {
            const auto& chunk = chunk_plan[static_cast<size_t>(idx)];
            const double local_cutoff = global_ub.load(std::memory_order_relaxed);
            double remaining_time_sec = std::numeric_limits<double>::infinity();
            if (std::isfinite(global_time_limit_sec)) {
                remaining_time_sec = global_time_limit_sec - paramip_timer.elapsed_seconds();
                if (remaining_time_sec <= 1e-6) {
                    SolveResult timeout;
                    timeout.status = SolveResult::Status::TimeLimit;
                    timeout.objective = 0.0;
                    timeout.bound = 0.0;
                    timeout.gap = 1.0;
                    timeout.time_seconds = 0.0;
                    chunk_results[static_cast<size_t>(idx)] = timeout;
                    return;
                }
            }
            const int32_t chunk_seed = selected_root_seed + 104729 * (chunk.id + 1);
            auto sub_opts = make_paramip_subsolve_options(
                options, encode_chunk_fixings(chunk), local_cutoff,
                remaining_time_sec, !output_flag,
                chunk_seed, std::nullopt);
            Model sub_model;
            sub_model.set_problem(problem_);
            auto r = sub_model.solve(sub_opts);
            chunk_results[static_cast<size_t>(idx)] = r;
            if (r.has_solution() && std::isfinite(r.objective)) {
                double prev = global_ub.load(std::memory_order_relaxed);
                while (r.objective + 1e-9 < prev) {
                    if (global_ub.compare_exchange_weak(
                            prev, r.objective,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        break;
                    }
                }
            }
        };

        if (deterministic) {
            for (int32_t i = 0; i < static_cast<int32_t>(chunk_plan.size()); ++i) {
                run_chunk(i);
            }
        } else {
            const int32_t workers = std::min<int32_t>(
                std::max<int32_t>(1, paramip_worker_target),
                static_cast<int32_t>(chunk_plan.size()));
            std::atomic<int32_t> next{0};
            std::vector<std::thread> pool;
            pool.reserve(static_cast<size_t>(workers));
            for (int32_t w = 0; w < workers; ++w) {
                pool.emplace_back([&] {
                    for (;;) {
                        const int32_t idx = next.fetch_add(1, std::memory_order_relaxed);
                        if (idx >= static_cast<int32_t>(chunk_plan.size())) break;
                        run_chunk(idx);
                    }
                });
            }
            for (auto& t : pool) t.join();
        }

        SolveResult out;
        out.status = SolveResult::Status::Error;
        out.objective = std::numeric_limits<double>::infinity();
        out.bound = std::numeric_limits<double>::infinity();
        out.gap = 1.0;
        out.time_seconds = paramip_timer.elapsed_seconds();

        bool any_error = false;
        bool any_timelimit = false;
        bool any_solution = false;
        bool all_infeasible = true;
        bool all_opt_or_infeasible = true;

        for (const auto& r : chunk_results) {
            out.nodes += r.nodes;
            out.total_cuts += r.total_cuts;
            out.separation_rounds += r.separation_rounds;
            out.dssr_epochs_enqueued += r.dssr_epochs_enqueued;
            out.dssr_epochs_committed += r.dssr_epochs_committed;
            out.dssr_checkpoint_count += r.dssr_checkpoint_count;
            out.dssr_commit_signature ^= r.dssr_commit_signature;
            aggregate_separator_stats(out.separator_stats, r.separator_stats);

            any_error = any_error || (r.status == SolveResult::Status::Error);
            any_timelimit = any_timelimit || (r.status == SolveResult::Status::TimeLimit);
            all_infeasible = all_infeasible && (r.status == SolveResult::Status::Infeasible);
            all_opt_or_infeasible = all_opt_or_infeasible &&
                (r.status == SolveResult::Status::Optimal ||
                 r.status == SolveResult::Status::Infeasible);
            if (std::isfinite(r.bound)) {
                out.bound = std::min(out.bound, r.bound);
            }
            if (r.has_solution() && std::isfinite(r.objective)) {
                any_solution = true;
                if (!std::isfinite(out.objective) || r.objective < out.objective - 1e-9) {
                    out.objective = r.objective;
                    out.tour = r.tour;
                    out.tour_arcs = r.tour_arcs;
                }
            }
        }

        if (any_solution) {
            if (!std::isfinite(out.bound)) out.bound = out.objective;
            if (std::abs(out.objective) > 1e-9) {
                out.gap = std::max(0.0, (out.objective - out.bound) / std::abs(out.objective));
            } else {
                out.gap = (out.bound <= out.objective + 1e-9) ? 0.0 : 1.0;
            }
            if (all_opt_or_infeasible) out.status = SolveResult::Status::Optimal;
            else if (any_timelimit) out.status = SolveResult::Status::TimeLimit;
            else out.status = SolveResult::Status::Feasible;
        } else if (all_infeasible) {
            out.status = SolveResult::Status::Infeasible;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 0.0;
        } else if (any_timelimit) {
            out.status = SolveResult::Status::TimeLimit;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 1.0;
        } else if (any_error) {
            out.status = SolveResult::Status::Error;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 1.0;
        } else {
            out.status = SolveResult::Status::Feasible;
        }

        logger_.log("ParaMIP static_root done: chunks={}, status={}, obj={}, bound={}, time={}s",
                    static_cast<int32_t>(chunk_plan.size()),
                    static_cast<int>(out.status),
                    out.objective,
                    out.bound,
                    out.time_seconds);
        return out;
    }

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

    // Suppress HiGHS console output; forward through our logger if enabled
    highs.setOptionValue("output_flag", false);
    if (output_flag) {
        highs.setLogCallback(
            [](HighsLogType /*type*/, const char* message, void* data) {
                static_cast<Logger*>(data)->log(std::string_view{message});
            },
            &logger_);
    }

    // Staged preprocessing:
    // Stage 1: parallel source/target 2-cycle bounds + fast warm-start
    // Stage 2: adaptive second warm-start + optional all-pairs + deeper s-t ng/DSSR
    // correction = profit(source) when s == t (tour: depot profit double-subtracted)
    double correction = problem_.is_tour() ? problem_.profit(source) : 0.0;

    const int32_t m = problem_.num_edges();

    std::vector<double> fwd_bounds, bwd_bounds;
    std::vector<double> all_pairs;
    heuristic::HeuristicResult warm_start{{}, std::numeric_limits<double>::infinity()};
    double warm_start_ub = std::numeric_limits<double>::infinity();
    int32_t ng_size_used = 1;
    double ng_path_ub = std::numeric_limits<double>::infinity();
    std::vector<int32_t> ng_path_nodes;
    int32_t stage1_elim_edges = 0;
    double stage1_elim_ratio = 0.0;

    // Use cutoff as UB when provided (known optimal for prove-only mode)
    if (cutoff < std::numeric_limits<double>::infinity()) {
        warm_start_ub = cutoff;
    }

    if (internal_skip_preproc && prepared_context) {
        fwd_bounds = prepared_context->fwd_bounds;
        bwd_bounds = prepared_context->bwd_bounds;
        if (all_pairs_propagation) {
            all_pairs = prepared_context->all_pairs;
        }
        warm_start = prepared_context->warm_start;
        warm_start_ub = prepared_context->warm_start_ub;
        if (cutoff < std::numeric_limits<double>::infinity()) {
            warm_start_ub = std::min(warm_start_ub, cutoff);
        }
        ng_size_used = prepared_context->ng_size_used;
        correction = prepared_context->correction;
        stage1_elim_edges = prepared_context->stage1_elim_edges;
        stage1_elim_ratio = prepared_context->stage1_elim_ratio;
        logger_.log(
            "Preprocessing: 0s, UB={}{} [stage1={}, edge-elim estimate: {}/{} ({:.1f}%), reused]",
            warm_start_ub,
            (cutoff < std::numeric_limits<double>::infinity() ? " (cutoff)" : ""),
            prepared_context->stage1_backend_name,
            stage1_elim_edges, m, 100.0 * stage1_elim_ratio);
    } else {
        Timer preproc_timer;

        // Stage 1: fast startup bounds + first incumbent
        heuristic::HeuristicResult fast_start{{}, std::numeric_limits<double>::infinity()};
        enum class Stage1BoundsBackend {
            two_cycle,
            ng1,
        };
        Stage1BoundsBackend stage1_backend = Stage1BoundsBackend::two_cycle;
        if (preproc_stage1_bounds == "ng1") {
            stage1_backend = Stage1BoundsBackend::ng1;
        } else if (preproc_stage1_bounds == "auto") {
            // Favor the lightweight 2-cycle kernel on medium/large instances.
            stage1_backend = (n <= 40) ? Stage1BoundsBackend::ng1
                                       : Stage1BoundsBackend::two_cycle;
        }
        const char* stage1_backend_name =
            (stage1_backend == Stage1BoundsBackend::ng1) ? "ng1" : "two_cycle";
        if (workflow_dump) {
            logger_.log("Workflow DAG (startup):");
            logger_.log("  stage1_bounds({}) + fast_warm_start -> ub0", stage1_backend_name);
            logger_.log("  ub0 + stage1_bounds -> stage1_edge_elim_ratio -> second_warm_start(adaptive)");
            logger_.log("  stage2 parallel: second_warm_start | ng_refinement | all_pairs(optional)");
        }

        tbb::task_group stage1_tg;
        if (stage1_backend == Stage1BoundsBackend::two_cycle) {
            stage1_tg.run([&] {
                fwd_bounds = preprocess::labeling_from(problem_, source);
            });
            if (source != target) {
                stage1_tg.run([&] {
                    bwd_bounds = preprocess::labeling_from(problem_, target);
                });
            }
            ng_size_used = 1;
        } else {
            auto fast_stage_opts = dssr_options;
            fast_stage_opts.initial_ng_size = 1;
            fast_stage_opts.max_ng_size = 1;
            fast_stage_opts.dssr_iterations = 1;
            stage1_tg.run([&] {
                auto fast_bounds = preprocess::ng::compute_bounds(
                    problem_, source, target, fast_stage_opts);
                fwd_bounds = std::move(fast_bounds.fwd);
                bwd_bounds = std::move(fast_bounds.bwd);
                ng_size_used = fast_bounds.ng_size;
                if (fast_bounds.elementary_path_found) {
                    ng_path_ub = fast_bounds.elementary_path_cost;
                    ng_path_nodes = std::move(fast_bounds.elementary_path);
                }
            });
        }
        if (!disable_heuristics) {
            const int fast_restarts = std::max<int32_t>(3, preproc_fast_restarts);
            const double fast_budget_ms = deterministic ? 0.0 : preproc_fast_budget_ms;
            stage1_tg.run([&] {
                fast_start = heuristic::build_initial_solution(
                    problem_, fast_restarts, fast_budget_ms, non_highs_work_budget);
            });
        }
        stage1_tg.wait();
        if (stage1_backend == Stage1BoundsBackend::two_cycle && source == target) {
            bwd_bounds = fwd_bounds;
        }
        if (!disable_heuristics) {
            warm_start = std::move(fast_start);
        }

        // Initial UB from cutoff and/or fast warm-start
        if (cutoff >= std::numeric_limits<double>::infinity()) {
            warm_start_ub = warm_start.objective;
        } else if (std::isfinite(warm_start.objective)) {
            warm_start_ub = std::min(warm_start_ub, warm_start.objective);
        }

        // First elimination estimate drives adaptive stage-2 warm-start policy
        if (edge_elimination
            && std::isfinite(warm_start_ub)
            && !fwd_bounds.empty() && !bwd_bounds.empty()) {
            auto stage1_eliminated = preprocess::edge_elimination(
                problem_, fwd_bounds, bwd_bounds, warm_start_ub, correction);
            stage1_elim_edges = static_cast<int32_t>(
                std::count(stage1_eliminated.begin(), stage1_eliminated.end(), true));
            stage1_elim_ratio = (m > 0)
                ? static_cast<double>(stage1_elim_edges) / static_cast<double>(m)
                : 0.0;
        }

        bool run_second_ws = !disable_heuristics;
        if (run_second_ws && preproc_adaptive) {
            const double min_elim = (n >= preproc_second_ws_large_n)
                ? preproc_second_ws_min_elim_large
                : preproc_second_ws_min_elim;
            run_second_ws = stage1_elim_ratio >= min_elim;
        }

        heuristic::HeuristicResult second_start{{}, std::numeric_limits<double>::infinity()};
        std::vector<double> refined_fwd;
        std::vector<double> refined_bwd;
        int32_t refined_ng_size = ng_size_used;
        bool refined_path_found = false;
        double refined_path_ub = std::numeric_limits<double>::infinity();
        std::vector<int32_t> refined_path_nodes;
        bool refined_bounds_ready = false;

        const bool run_ng_refinement =
            (dssr_options.dssr_iterations > 1
             || dssr_options.max_ng_size > 1
             || dssr_options.initial_ng_size > 1);

        tbb::task_group stage2_tg;
        if (run_second_ws) {
            int second_restarts = std::max(preproc_fast_restarts, std::clamp(n, 20, 200));
            double second_budget_ms = 0.0;
            if (!deterministic) {
                const double scaled = preproc_second_ws_budget_ms_min
                    + preproc_second_ws_budget_scale * static_cast<double>(n) * stage1_elim_ratio;
                second_budget_ms = std::clamp(
                    scaled, preproc_second_ws_budget_ms_min, preproc_second_ws_budget_ms_max);
            }
            stage2_tg.run([&] {
                second_start = heuristic::build_initial_solution(
                    problem_, second_restarts, second_budget_ms, non_highs_work_budget,
                    fwd_bounds, bwd_bounds, correction, warm_start_ub);
            });
        }

        if (all_pairs_propagation) {
            stage2_tg.run([&] {
                constexpr double inf = std::numeric_limits<double>::infinity();
                all_pairs.assign(static_cast<size_t>(n) * n, inf);
                tbb::parallel_for(0, n, [&](int32_t s) {
                    auto row = preprocess::forward_labeling(problem_, s);
                    std::copy(row.begin(), row.end(),
                              all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
                });
            });
        }

        if (run_ng_refinement) {
            stage2_tg.run([&] {
                auto refined = preprocess::ng::compute_bounds(
                    problem_, source, target, dssr_options);
                refined_fwd = std::move(refined.fwd);
                refined_bwd = std::move(refined.bwd);
                refined_ng_size = refined.ng_size;
                refined_bounds_ready = true;
                refined_path_found = refined.elementary_path_found;
                refined_path_ub = refined.elementary_path_cost;
                refined_path_nodes = std::move(refined.elementary_path);
            });
        }
        stage2_tg.wait();

        if (!second_start.col_values.empty()
            && second_start.objective + 1e-9 < warm_start.objective) {
            warm_start = std::move(second_start);
        }

        if (refined_bounds_ready
            && refined_fwd.size() == static_cast<size_t>(n)
            && refined_bwd.size() == static_cast<size_t>(n)) {
            for (int32_t i = 0; i < n; ++i) {
                fwd_bounds[i] = std::max(fwd_bounds[i], refined_fwd[i]);
                bwd_bounds[i] = std::max(bwd_bounds[i], refined_bwd[i]);
            }
            ng_size_used = std::max(ng_size_used, refined_ng_size);
            if (refined_path_found && refined_path_ub + 1e-9 < ng_path_ub) {
                ng_path_ub = refined_path_ub;
                ng_path_nodes = std::move(refined_path_nodes);
            }
        } else if (all_pairs_propagation && !all_pairs.empty()) {
            // If ng refinement is disabled, use source/target all-pairs rows.
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

        // Heuristic UB only used when no cutoff was provided
        if (cutoff >= std::numeric_limits<double>::infinity()) {
            warm_start_ub = warm_start.objective;
        } else if (std::isfinite(warm_start.objective)) {
            warm_start_ub = std::min(warm_start_ub, warm_start.objective);
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
        logger_.log("Preprocessing: {}s, UB={}{} [stage1={}, edge-elim estimate: {}/{} ({:.1f}%)]",
                    preproc_timer.elapsed_seconds(),
                    warm_start_ub,
                    (cutoff < std::numeric_limits<double>::infinity() ? " (cutoff)" : ""),
                    stage1_backend_name,
                    stage1_elim_edges, m, 100.0 * stage1_elim_ratio);
    }

    if (paramip_mode == "static_root" && !chunk_plan.empty() && !internal_skip_preproc) {
        auto prepared_ctx = std::make_shared<PreparedSolveContext>();
        prepared_ctx->fwd_bounds = fwd_bounds;
        prepared_ctx->bwd_bounds = bwd_bounds;
        prepared_ctx->all_pairs = all_pairs;
        prepared_ctx->warm_start = warm_start;
        prepared_ctx->warm_start_ub = warm_start_ub;
        prepared_ctx->ng_size_used = ng_size_used;
        prepared_ctx->correction = correction;
        prepared_ctx->stage1_elim_edges = stage1_elim_edges;
        prepared_ctx->stage1_elim_ratio = stage1_elim_ratio;
        prepared_ctx->stage1_backend_name = "prepared";
        const int64_t prepared_ctx_id = register_prepared_context(prepared_ctx);
        struct PreparedCtxGuard {
            int64_t id = 0;
            ~PreparedCtxGuard() { unregister_prepared_context(id); }
        } prepared_ctx_guard{prepared_ctx_id};

        Timer paramip_timer;
        std::vector<SolveResult> chunk_results(chunk_plan.size());
        std::atomic<double> global_ub{warm_start_ub};
        double global_time_limit_sec = std::numeric_limits<double>::infinity();
        for (const auto& [k, v] : options) {
            if (k != "time_limit") continue;
            try {
                global_time_limit_sec = std::stod(v);
            } catch (...) {
                global_time_limit_sec = std::numeric_limits<double>::infinity();
            }
            break;
        }

        int32_t selected_root_seed = base_random_seed;
        const std::string root_pick_policy =
            (paramip_root_pick == "auto")
                ? (deterministic ? "best" : "first")
                : paramip_root_pick;
        const int32_t root_probe_count =
            (paramip_root_probes > 0)
                ? paramip_root_probes
                : std::min<int32_t>(std::max<int32_t>(2, paramip_worker_target), 8);

        if (root_probe_count > 1) {
            struct RootProbeResult {
                int32_t probe_idx = -1;
                int32_t seed = 0;
                int32_t completion_rank = std::numeric_limits<int32_t>::max();
                SolveResult result;
            };

            std::vector<RootProbeResult> probes(static_cast<size_t>(root_probe_count));
            for (int32_t i = 0; i < root_probe_count; ++i) {
                probes[static_cast<size_t>(i)].probe_idx = i;
                probes[static_cast<size_t>(i)].seed = base_random_seed + i;
            }

            double stage0_remaining = std::numeric_limits<double>::infinity();
            if (std::isfinite(global_time_limit_sec)) {
                stage0_remaining = global_time_limit_sec - paramip_timer.elapsed_seconds();
            }
            if (stage0_remaining > 1e-6 || !std::isfinite(stage0_remaining)) {
                const double stage0_probe_limit =
                    std::isfinite(stage0_remaining)
                        ? std::max(0.05, std::min(5.0, 0.25 * stage0_remaining))
                        : std::numeric_limits<double>::infinity();
                const std::string root_fixings = encode_paramip_fixings(fixed_y);
                std::atomic<int32_t> completion_counter{0};
                auto run_probe = [&](int32_t idx) {
                    auto& probe = probes[static_cast<size_t>(idx)];
                    auto sub_opts = make_paramip_subsolve_options(
                        options, root_fixings, global_ub.load(std::memory_order_relaxed),
                        stage0_probe_limit, true,
                        probe.seed, 1);
                    sub_opts.push_back({"_internal_skip_preproc", "true"});
                    sub_opts.push_back({"_internal_prepared_ctx_id",
                                        std::to_string(prepared_ctx_id)});
                    sub_opts.push_back({"_internal_disable_async_dssr", "true"});
                    Model sub_model;
                    sub_model.set_problem(problem_);
                    probe.result = sub_model.solve(sub_opts);
                    probe.completion_rank = completion_counter.fetch_add(
                        1, std::memory_order_relaxed);
                };

                const int32_t workers = std::min<int32_t>(
                    paramip_worker_target, root_probe_count);
                if (workers <= 1) {
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        run_probe(i);
                    }
                } else {
                    std::atomic<int32_t> next{0};
                    std::vector<std::thread> pool;
                    pool.reserve(static_cast<size_t>(workers));
                    for (int32_t w = 0; w < workers; ++w) {
                        pool.emplace_back([&] {
                            for (;;) {
                                const int32_t idx =
                                    next.fetch_add(1, std::memory_order_relaxed);
                                if (idx >= root_probe_count) break;
                                run_probe(idx);
                            }
                        });
                    }
                    for (auto& t : pool) t.join();
                }

                auto has_finite_bound = [](const RootProbeResult& p) {
                    return std::isfinite(p.result.bound);
                };
                auto has_incumbent = [](const RootProbeResult& p) {
                    if (!std::isfinite(p.result.objective)) return false;
                    if (p.result.status == SolveResult::Status::Optimal ||
                        p.result.status == SolveResult::Status::Feasible) {
                        return true;
                    }
                    if (p.result.status == SolveResult::Status::TimeLimit) {
                        return !p.result.tour.empty() || !p.result.tour_arcs.empty();
                    }
                    return false;
                };

                int32_t chosen = -1;
                if (root_pick_policy == "first") {
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        const auto& p = probes[static_cast<size_t>(i)];
                        if (p.result.status == SolveResult::Status::Error) continue;
                        if (!has_finite_bound(p)) continue;
                        if (chosen < 0 ||
                            p.completion_rank <
                                probes[static_cast<size_t>(chosen)].completion_rank) {
                            chosen = i;
                        }
                    }
                    if (chosen < 0) {
                        for (int32_t i = 0; i < root_probe_count; ++i) {
                            const auto& p = probes[static_cast<size_t>(i)];
                            if (p.result.status == SolveResult::Status::Error) continue;
                            if (chosen < 0 ||
                                p.completion_rank <
                                    probes[static_cast<size_t>(chosen)].completion_rank) {
                                chosen = i;
                            }
                        }
                    }
                } else {
                    for (int32_t i = 0; i < root_probe_count; ++i) {
                        const auto& p = probes[static_cast<size_t>(i)];
                        if (p.result.status == SolveResult::Status::Error) continue;
                        if (chosen < 0) {
                            chosen = i;
                            continue;
                        }
                        const auto& best = probes[static_cast<size_t>(chosen)];
                        const bool p_has_bound = has_finite_bound(p);
                        const bool b_has_bound = has_finite_bound(best);
                        if (p_has_bound && !b_has_bound) {
                            chosen = i;
                            continue;
                        }
                        if (p_has_bound && b_has_bound &&
                            p.result.bound > best.result.bound + 1e-9) {
                            chosen = i;
                            continue;
                        }
                        if (p_has_bound && b_has_bound &&
                            std::abs(p.result.bound - best.result.bound) <= 1e-9) {
                            const bool p_has_sol = has_incumbent(p);
                            const bool b_has_sol = has_incumbent(best);
                            if (p_has_sol && !b_has_sol) {
                                chosen = i;
                                continue;
                            }
                            if (p_has_sol && b_has_sol &&
                                p.result.objective < best.result.objective - 1e-9) {
                                chosen = i;
                                continue;
                            }
                        }
                    }
                }

                if (chosen >= 0) {
                    const auto& picked = probes[static_cast<size_t>(chosen)];
                    selected_root_seed = picked.seed;
                    double best_probe_obj = std::numeric_limits<double>::infinity();
                    for (const auto& p : probes) {
                        if (!has_incumbent(p)) continue;
                        best_probe_obj = std::min(best_probe_obj, p.result.objective);
                    }
                    if (std::isfinite(best_probe_obj)) {
                        double prev = global_ub.load(std::memory_order_relaxed);
                        while (best_probe_obj + 1e-9 < prev) {
                            if (global_ub.compare_exchange_weak(
                                    prev, best_probe_obj,
                                    std::memory_order_relaxed,
                                    std::memory_order_relaxed)) {
                                break;
                            }
                        }
                    }
                    logger_.log(
                        "ParaMIP stage0: probes={}, policy={}, chosen_probe={}, seed={}, bound={}, obj={}, status={}",
                        root_probe_count, root_pick_policy, chosen, selected_root_seed,
                        picked.result.bound, picked.result.objective,
                        static_cast<int>(picked.result.status));
                } else {
                    logger_.log(
                        "ParaMIP stage0: probes={}, policy={}, no valid root candidate (keeping seed={})",
                        root_probe_count, root_pick_policy, selected_root_seed);
                }
            }
        }

        auto encode_chunk_fixings = [&](const ParamipChunkPlan& chunk) {
            std::vector<int8_t> merged = fixed_y;
            if (merged.empty()) {
                merged.assign(static_cast<size_t>(n), static_cast<int8_t>(-1));
            }
            for (const auto& f : chunk.fixings) {
                if (f.node < 0 || f.node >= n) continue;
                merged[static_cast<size_t>(f.node)] = static_cast<int8_t>(f.value ? 1 : 0);
            }
            return encode_paramip_fixings(merged);
        };

        auto run_chunk = [&](int32_t idx) {
            const auto& chunk = chunk_plan[static_cast<size_t>(idx)];
            const double local_cutoff = global_ub.load(std::memory_order_relaxed);
            double remaining_time_sec = std::numeric_limits<double>::infinity();
            if (std::isfinite(global_time_limit_sec)) {
                remaining_time_sec = global_time_limit_sec - paramip_timer.elapsed_seconds();
                if (remaining_time_sec <= 1e-6) {
                    SolveResult timeout;
                    timeout.status = SolveResult::Status::TimeLimit;
                    timeout.objective = 0.0;
                    timeout.bound = 0.0;
                    timeout.gap = 1.0;
                    timeout.time_seconds = 0.0;
                    chunk_results[static_cast<size_t>(idx)] = timeout;
                    return;
                }
            }
            const int32_t chunk_seed = selected_root_seed + 104729 * (chunk.id + 1);
            auto sub_opts = make_paramip_subsolve_options(
                options, encode_chunk_fixings(chunk), local_cutoff,
                remaining_time_sec, !output_flag,
                chunk_seed, std::nullopt);
            sub_opts.push_back({"_internal_skip_preproc", "true"});
            sub_opts.push_back({"_internal_prepared_ctx_id",
                                std::to_string(prepared_ctx_id)});
            Model sub_model;
            sub_model.set_problem(problem_);
            auto r = sub_model.solve(sub_opts);
            chunk_results[static_cast<size_t>(idx)] = r;
            if (r.has_solution() && std::isfinite(r.objective)) {
                double prev = global_ub.load(std::memory_order_relaxed);
                while (r.objective + 1e-9 < prev) {
                    if (global_ub.compare_exchange_weak(
                            prev, r.objective,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        break;
                    }
                }
            }
        };

        if (deterministic) {
            for (int32_t i = 0; i < static_cast<int32_t>(chunk_plan.size()); ++i) {
                run_chunk(i);
            }
        } else {
            const int32_t workers = std::min<int32_t>(
                std::max<int32_t>(1, paramip_worker_target),
                static_cast<int32_t>(chunk_plan.size()));
            std::atomic<int32_t> next{0};
            std::vector<std::thread> pool;
            pool.reserve(static_cast<size_t>(workers));
            for (int32_t w = 0; w < workers; ++w) {
                pool.emplace_back([&] {
                    for (;;) {
                        const int32_t idx = next.fetch_add(1, std::memory_order_relaxed);
                        if (idx >= static_cast<int32_t>(chunk_plan.size())) break;
                        run_chunk(idx);
                    }
                });
            }
            for (auto& t : pool) t.join();
        }

        SolveResult out;
        out.status = SolveResult::Status::Error;
        out.objective = std::numeric_limits<double>::infinity();
        out.bound = std::numeric_limits<double>::infinity();
        out.gap = 1.0;
        out.time_seconds = paramip_timer.elapsed_seconds();

        bool any_error = false;
        bool any_timelimit = false;
        bool any_solution = false;
        bool all_infeasible = true;
        bool all_opt_or_infeasible = true;

        for (const auto& r : chunk_results) {
            out.nodes += r.nodes;
            out.total_cuts += r.total_cuts;
            out.separation_rounds += r.separation_rounds;
            out.dssr_epochs_enqueued += r.dssr_epochs_enqueued;
            out.dssr_epochs_committed += r.dssr_epochs_committed;
            out.dssr_checkpoint_count += r.dssr_checkpoint_count;
            out.dssr_commit_signature ^= r.dssr_commit_signature;
            aggregate_separator_stats(out.separator_stats, r.separator_stats);

            any_error = any_error || (r.status == SolveResult::Status::Error);
            any_timelimit = any_timelimit || (r.status == SolveResult::Status::TimeLimit);
            all_infeasible = all_infeasible && (r.status == SolveResult::Status::Infeasible);
            all_opt_or_infeasible = all_opt_or_infeasible &&
                (r.status == SolveResult::Status::Optimal ||
                 r.status == SolveResult::Status::Infeasible);
            if (std::isfinite(r.bound)) {
                out.bound = std::min(out.bound, r.bound);
            }
            if (r.has_solution() && std::isfinite(r.objective)) {
                any_solution = true;
                if (!std::isfinite(out.objective) || r.objective < out.objective - 1e-9) {
                    out.objective = r.objective;
                    out.tour = r.tour;
                    out.tour_arcs = r.tour_arcs;
                }
            }
        }

        if (any_solution) {
            if (!std::isfinite(out.bound)) out.bound = out.objective;
            if (std::abs(out.objective) > 1e-9) {
                out.gap = std::max(0.0, (out.objective - out.bound) / std::abs(out.objective));
            } else {
                out.gap = (out.bound <= out.objective + 1e-9) ? 0.0 : 1.0;
            }
            if (all_opt_or_infeasible) out.status = SolveResult::Status::Optimal;
            else if (any_timelimit) out.status = SolveResult::Status::TimeLimit;
            else out.status = SolveResult::Status::Feasible;
        } else if (all_infeasible) {
            out.status = SolveResult::Status::Infeasible;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 0.0;
        } else if (any_timelimit) {
            out.status = SolveResult::Status::TimeLimit;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 1.0;
        } else if (any_error) {
            out.status = SolveResult::Status::Error;
            out.objective = 0.0;
            if (!std::isfinite(out.bound)) out.bound = 0.0;
            out.gap = 1.0;
        } else {
            out.status = SolveResult::Status::Feasible;
        }

        logger_.log("ParaMIP static_root done: chunks={}, status={}, obj={}, bound={}, time={}s",
                    static_cast<int32_t>(chunk_plan.size()),
                    static_cast<int>(out.status),
                    out.objective,
                    out.bound,
                    out.time_seconds);
        return out;
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
    bridge.set_heuristic_deterministic_mode(deterministic);
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

    const int32_t async_start_ng =
        std::max(ng_size_used + 1, dssr_options.initial_ng_size + 1);
    int32_t async_end_ng = dssr_options.max_ng_size;
    if (dssr_background_max_epochs > 0) {
        const int64_t capped_end =
            static_cast<int64_t>(async_start_ng) +
            static_cast<int64_t>(dssr_background_max_epochs) - 1;
        async_end_ng = std::min<int32_t>(
            async_end_ng,
            static_cast<int32_t>(
                std::min<int64_t>(capped_end, std::numeric_limits<int32_t>::max())));
    }
    const bool dssr_background_updates_enabled =
        dssr_background_updates &&
        !internal_disable_async_dssr &&
        !all_pairs_propagation &&
        async_end_ng >= async_start_ng;
    if (workflow_dump) {
        logger_.log("Workflow DAG (solve):");
        logger_.log("  build_formulation -> highs.run() with callbacks");
        if (dssr_background_updates_enabled) {
            logger_.log("  highs.run() || async_dssr_updates (ng_size {}..{})",
                        async_start_ng, async_end_ng);
        } else {
            logger_.log("  async_dssr_updates disabled for this run");
        }
        logger_.log("  paramip-style worker DAG: enabled via paramip_mode=static_root");
    }
    const bool deterministic_async_strict =
        dssr_background_updates_enabled && deterministic;
    const bool dssr_auto_policy =
        dssr_background_updates_enabled && dssr_background_policy == "auto";
    auto dssr_epoch_queue = deterministic_async_strict
        ? std::make_shared<model_detail::DssrEpochQueue>(async_start_ng)
        : std::shared_ptr<model_detail::DssrEpochQueue>{};
    auto dssr_checkpoint_clock = std::make_shared<model_detail::DeterministicCheckpointClock>();
    auto dssr_commit_mutex = std::make_shared<std::mutex>();
    auto dssr_request_mutex = std::make_shared<std::mutex>();
    auto dssr_request_cv = std::make_shared<std::condition_variable>();
    auto dssr_requested_epochs = std::make_shared<int32_t>(0);
    std::atomic<int64_t> dssr_epochs_enqueued{0};
    int64_t dssr_epochs_committed = 0;
    uint64_t dssr_commit_signature = kDssrSigOffset;

    auto commit_dssr_epoch = [&](model_detail::DssrEpochUpdate update) {
        bool bounds_changed =
            shared_bounds->publish_tightening(update.fwd, update.bwd, update.ng_size);
        bool ub_changed = false;
        bool incumbent_changed = false;

        if (update.elementary_path_found) {
            double prev = async_upper_bound->load(std::memory_order_relaxed);
            while (update.elementary_path_cost + 1e-9 < prev) {
                if (async_upper_bound->compare_exchange_weak(
                        prev, update.elementary_path_cost,
                        std::memory_order_relaxed,
                        std::memory_order_relaxed)) {
                    ub_changed = true;
                    break;
                }
            }
            if (!update.elementary_path.empty()) {
                auto start = build_ng_path_start(problem_, update.elementary_path);
                if (!start.col_values.empty()) {
                    incumbent_changed = async_incumbent_store->publish_if_better(
                        std::move(start.col_values), start.objective);
                }
            }
        }

        dssr_epochs_committed++;
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, static_cast<uint64_t>(update.epoch));
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, static_cast<uint64_t>(update.ng_size));
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, static_cast<uint64_t>(bounds_changed ? 1 : 0));
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, static_cast<uint64_t>(ub_changed ? 1 : 0));
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, static_cast<uint64_t>(incumbent_changed ? 1 : 0));
        dssr_commit_signature = dssr_sig_mix(
            dssr_commit_signature, dssr_sig_quantize(update.elementary_path_cost));
    };

    if (deterministic_async_strict && dssr_epoch_queue) {
        bridge.set_deterministic_checkpoint_hook(
            [dssr_epoch_queue, dssr_checkpoint_clock, dssr_commit_mutex,
             dssr_request_mutex, dssr_request_cv, dssr_requested_epochs,
             &commit_dssr_epoch, deterministic_async_strict]() {
                {
                    std::lock_guard<std::mutex> req_lock(*dssr_request_mutex);
                    (*dssr_requested_epochs)++;
                }
                dssr_request_cv->notify_one();
                std::lock_guard<std::mutex> lock(*dssr_commit_mutex);
                dssr_checkpoint_clock->checkpoint();
                if (deterministic_async_strict) {
                    dssr_epoch_queue->commit_next_blocking(
                        [&commit_dssr_epoch](model_detail::DssrEpochUpdate update) {
                            commit_dssr_epoch(std::move(update));
                        });
                } else {
                    dssr_epoch_queue->commit_ready(
                        [&commit_dssr_epoch](model_detail::DssrEpochUpdate update) {
                            commit_dssr_epoch(std::move(update));
                        });
                }
            });
    } else if (dssr_auto_policy) {
        bridge.set_deterministic_checkpoint_hook(
            [dssr_checkpoint_clock]() { dssr_checkpoint_clock->checkpoint(); });
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
    if (dssr_background_updates_enabled) {
        async_dssr_worker = std::thread([&]() {
            auto auto_snapshot = shared_bounds->snapshot();
            std::vector<double> auto_best_fwd = auto_snapshot.fwd;
            std::vector<double> auto_best_bwd = auto_snapshot.bwd;
            double auto_best_ub = async_upper_bound->load(std::memory_order_relaxed);
            model_detail::DssrAutoStopTracker auto_stop_tracker(
                dssr_background_auto_min_epochs,
                dssr_background_auto_no_progress_limit);
            auto mark_bounds_tightening = [&](std::vector<double>& current,
                                              const std::vector<double>& next) {
                if (current.size() != next.size()) return false;
                bool changed = false;
                for (size_t i = 0; i < current.size(); ++i) {
                    const double cur = current[i];
                    const double nxt = next[i];
                    if (!std::isfinite(nxt)) continue;
                    if (!std::isfinite(cur) || nxt > cur + 1e-12) {
                        current[i] = nxt;
                        changed = true;
                    }
                }
                return changed;
            };
            int32_t produced_epochs = 0;
            for (int32_t ng_size = async_start_ng;
                 ng_size <= async_end_ng;
                 ++ng_size) {
                if (deterministic_async_strict) {
                    std::unique_lock<std::mutex> req_lock(*dssr_request_mutex);
                    dssr_request_cv->wait(req_lock, [&]() {
                        return stop_async_dssr.load(std::memory_order_relaxed) ||
                               produced_epochs < *dssr_requested_epochs;
                    });
                    if (stop_async_dssr.load(std::memory_order_relaxed)) {
                        break;
                    }
                }
                if (!deterministic_async_strict &&
                    stop_async_dssr.load(std::memory_order_relaxed)) {
                    break;
                }
                if (non_highs_work_budget &&
                    !non_highs_work_budget->try_consume(1)) {
                    break;
                }
                auto stage_opts = dssr_options;
                stage_opts.initial_ng_size = ng_size;
                stage_opts.max_ng_size = ng_size;
                stage_opts.dssr_iterations = std::max<int32_t>(1, dssr_options.dssr_iterations / 2);
                auto bounds = preprocess::ng::compute_bounds(problem_, source, target, stage_opts);
                model_detail::DssrEpochUpdate update;
                update.epoch = ng_size;
                update.ng_size = bounds.ng_size;
                update.fwd = std::move(bounds.fwd);
                update.bwd = std::move(bounds.bwd);
                update.elementary_path_found = bounds.elementary_path_found;
                update.elementary_path_cost = bounds.elementary_path_cost;
                update.elementary_path = std::move(bounds.elementary_path);
                bool bounds_tightened_for_auto = false;
                bool ub_tightened_for_auto = false;
                if (dssr_auto_policy) {
                    bounds_tightened_for_auto |=
                        mark_bounds_tightening(auto_best_fwd, update.fwd);
                    bounds_tightened_for_auto |=
                        mark_bounds_tightening(auto_best_bwd, update.bwd);
                    if (update.elementary_path_found &&
                        update.elementary_path_cost + 1e-9 < auto_best_ub) {
                        auto_best_ub = update.elementary_path_cost;
                        ub_tightened_for_auto = true;
                    }
                }
                bool stage_improved = false;
                if (deterministic_async_strict) {
                    dssr_epoch_queue->enqueue(std::move(update));
                    dssr_epochs_enqueued.fetch_add(1, std::memory_order_relaxed);
                    if (dssr_auto_policy) {
                        stage_improved =
                            bounds_tightened_for_auto || ub_tightened_for_auto;
                    }
                } else {
                    bool bounds_changed =
                        shared_bounds->publish_tightening(update.fwd, update.bwd, update.ng_size);
                    bool ub_changed = false;
                    bool incumbent_changed = false;
                    if (update.elementary_path_found) {
                        double prev = async_upper_bound->load(std::memory_order_relaxed);
                        while (update.elementary_path_cost + 1e-9 < prev) {
                            if (async_upper_bound->compare_exchange_weak(
                                    prev, update.elementary_path_cost,
                                    std::memory_order_relaxed,
                                    std::memory_order_relaxed)) {
                                ub_changed = true;
                                break;
                            }
                        }
                        if (!update.elementary_path.empty()) {
                            auto start = build_ng_path_start(problem_, update.elementary_path);
                            if (!start.col_values.empty()) {
                                incumbent_changed = async_incumbent_store->publish_if_better(
                                    std::move(start.col_values), start.objective);
                            }
                        }
                    }
                    if (dssr_auto_policy) {
                        stage_improved = bounds_tightened_for_auto ||
                                         ub_changed || incumbent_changed;
                    } else {
                        stage_improved = bounds_changed || ub_changed || incumbent_changed;
                    }
                }

                produced_epochs++;
                if (dssr_auto_policy) {
                    const bool checkpoint_active =
                        dssr_checkpoint_clock->value() > 0;
                    if (auto_stop_tracker.observe_stage(
                            stage_improved, checkpoint_active, produced_epochs)) {
                        break;
                    }
                }
            }
            if (deterministic_async_strict) {
                dssr_epoch_queue->mark_producer_done();
            }
        });
    }

    highs.run();

    stop_async_dssr.store(true, std::memory_order_relaxed);
    if (deterministic_async_strict) {
        dssr_request_cv->notify_all();
    }
    if (async_dssr_worker.joinable()) {
        async_dssr_worker.join();
    }
    if (deterministic_async_strict && dssr_epoch_queue) {
        std::lock_guard<std::mutex> lock(*dssr_commit_mutex);
        dssr_checkpoint_clock->checkpoint();
        dssr_epoch_queue->commit_all_ordered(
            [&commit_dssr_epoch](model_detail::DssrEpochUpdate update) {
                commit_dssr_epoch(std::move(update));
            });
    }

    // Clear static callbacks/state to avoid leaking between solves
    HighsUserSeparator::clearCallback();

    auto result = bridge.extract_result();
    result.time_seconds = timer.elapsed_seconds();
    result.dssr_epochs_enqueued =
        dssr_epochs_enqueued.load(std::memory_order_relaxed);
    result.dssr_epochs_committed = dssr_epochs_committed;
    result.dssr_checkpoint_count = dssr_checkpoint_clock->value();
    result.dssr_commit_signature =
        (dssr_epochs_committed > 0) ? dssr_commit_signature : 0;
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
