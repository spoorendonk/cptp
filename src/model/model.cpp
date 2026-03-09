#include "model/model.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string_view>
#include <thread>

#include "model/highs_bridge.h"
#include "parallel/parallel.h"

#ifdef __linux__
#include <unistd.h>
#endif

#include "heuristic/primal_heuristic.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsUserSeparator.h"
#include "preprocess/edge_elimination.h"
#include "sep/comb_separator.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/rglm_separator.h"
#include "sep/sec_separator.h"
#include "sep/spi_separator.h"
#include "util/timer.h"

namespace cptp {

namespace {

int32_t count_finite_bounds(std::span<const double> bounds) {
  return static_cast<int32_t>(std::count_if(
      bounds.begin(), bounds.end(), [](double v) { return std::isfinite(v); }));
}

bool has_feasible_heuristic(const heuristic::HeuristicResult& r) {
  return !r.col_values.empty() && std::isfinite(r.objective) &&
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

constexpr int32_t kWarmStartLsStartsPerThread = 4;

void configure_hyperplane_branching(const Problem& problem, HiGHSBridge& bridge,
                                    Logger& logger, std::string& branch_hyper,
                                    int32_t hyper_sb_max_depth,
                                    int32_t hyper_sb_iter_limit,
                                    int32_t hyper_sb_min_reliable,
                                    int32_t hyper_sb_max_candidates) {
  // Hyperplane branching: dynamic constraint branching
  static constexpr std::string_view valid_modes[] = {
      "off", "pairs", "clusters", "demand", "cardinality", "all"};
  if (std::find(std::begin(valid_modes), std::end(valid_modes), branch_hyper) ==
      std::end(valid_modes)) {
    logger.log("Warning: unknown branch_hyper mode '{}', using 'off'",
               branch_hyper);
    branch_hyper = "off";
  }
  if (branch_hyper == "off") {
    // Ensure no stale global hyperplane callback leaks from prior solves.
    HighsUserSeparator::setBranchingCallback({});
    HighsUserSeparator::clearStack();
    return;
  }

  HighsUserSeparator::StrongBranchConfig sb_cfg;
  sb_cfg.max_depth = hyper_sb_max_depth;
  sb_cfg.iter_limit = hyper_sb_iter_limit;
  sb_cfg.min_reliable = hyper_sb_min_reliable;
  sb_cfg.max_sb_candidates = hyper_sb_max_candidates;
  HighsUserSeparator::setStrongBranchConfig(sb_cfg);

  const int32_t n = bridge.num_nodes();
  const int32_t y_off = bridge.y_offset();
  const auto& graph = problem.graph();

  std::vector<double> demands(n, 0.0);
  for (int32_t i = 0; i < n; ++i) demands[i] = problem.demand(i);

  // Ryan-Foster pairs: for each non-depot node, pair with nearest neighbor.
  struct Pair {
    int32_t i, j;
  };
  std::vector<Pair> rf_pairs;
  if (branch_hyper == "pairs" || branch_hyper == "all") {
    std::vector<bool> paired(n, false);
    for (int32_t i = 0; i < n; ++i) {
      if (i == problem.source()) continue;
      double best_cost = std::numeric_limits<double>::infinity();
      int32_t best_j = -1;
      for (auto e : graph.incident_edges(i)) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        int32_t nb = (u == i) ? v : u;
        if (nb == problem.source()) continue;
        double c = problem.edge_cost(e);
        if (c < best_cost) {
          best_cost = c;
          best_j = nb;
        }
      }
      if (best_j >= 0) {
        // Deduplicate: only add (min, max) pairs.
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

  // Greedy nearest-neighbor clusters of size ~4.
  std::vector<std::vector<int32_t>> clusters;
  if (branch_hyper == "clusters" || branch_hyper == "all") {
    constexpr int32_t CLUSTER_SIZE = 4;
    std::vector<bool> assigned(n, false);
    assigned[problem.source()] = true;

    // Sort nodes by demand (descending) to seed clusters.
    std::vector<int32_t> nodes_by_demand;
    for (int32_t i = 0; i < n; ++i) {
      if (i != problem.source()) nodes_by_demand.push_back(i);
    }
    std::sort(nodes_by_demand.begin(), nodes_by_demand.end(),
              [&](int32_t a, int32_t b) { return demands[a] > demands[b]; });

    for (int32_t seed : nodes_by_demand) {
      if (assigned[seed]) continue;
      std::vector<int32_t> cluster = {seed};
      assigned[seed] = true;

      // Add nearest unassigned neighbors.
      while (static_cast<int32_t>(cluster.size()) < CLUSTER_SIZE) {
        double best_cost = std::numeric_limits<double>::infinity();
        int32_t best_nb = -1;
        for (int32_t ci : cluster) {
          for (auto e : graph.incident_edges(ci)) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            int32_t nb = (u == ci) ? v : u;
            if (assigned[nb]) continue;
            double c = problem.edge_cost(e);
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
      [y_off, n, source = problem.source(), demands = std::move(demands),
       branch_hyper, rf_pairs = std::move(rf_pairs),
       clusters = std::move(clusters)](const HighsMipSolver& mipsolver)
          -> std::vector<HighsUserSeparator::HyperplaneCandidate> {
        const auto& sol = mipsolver.mipdata_->lp.getSolution().col_value;
        std::vector<HighsUserSeparator::HyperplaneCandidate> result;

        if (branch_hyper == "pairs" || branch_hyper == "all") {
          for (auto& [pi, pj] : rf_pairs) {
            double yi = sol[y_off + pi], yj = sol[y_off + pj];
            // Skip if either variable is near-integer.
            if (yi < 0.05 || yi > 0.95 || yj < 0.05 || yj > 0.95) continue;
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
            // Cluster demand hyperplane: sum(demand_i * y_i) for i in cluster.
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
              if (frac > 0.01 && frac < 0.99) result.push_back(std::move(c));
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

}  // namespace

Model::Model() = default;

void Model::set_problem(Problem prob) {
  problem_ = std::move(prob);
  built_ = true;
}

void Model::set_graph(int32_t num_nodes, std::span<const Edge> edges,
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

void Model::add_capacity_resource(std::span<const double> demands,
                                  double limit) {
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
                 std::move(profits_), std::move(demands), capacity, source_,
                 target_);
  built_ = true;
}

SolveResult Model::solve(const SolverOptions& options) {
  build_problem();

  Timer timer;

  int32_t random_seed = 0;
  bool output_flag = true;
  for (const auto& [key, value] : options) {
    if (key == "output_flag") {
      output_flag = (value == "true" || value == "1");
    }
  }
  logger_.set_enabled(output_flag);
  logger_.log("cptp v0.1");

  Highs highs;
  // Our defaults (user options can override)
  highs.setOptionValue("presolve", "off");

  // Intercept our custom options, forward the rest to HiGHS
  // --- Bounds & Preprocessing ---
  bool labeling_elim_bounds = true;
  bool all_pairs_bounds = false;
  bool edge_elimination = true;
  bool edge_elimination_nodes = true;
  bool bounds_propagation = true;

  // --- Separation (global) ---
  int32_t separation_interval = 1;
  double separation_tol = sep::kDefaultFracTol;
  int32_t max_cuts_per_round = 10;

  // --- Heuristics: warm-start ---
  bool heu_ws = true;
  int32_t heu_ws_ls_max_iter = 1000;

  // --- Heuristics: LP-guided ---
  bool heu_lpg = true;
  int heu_lpg_strategy = 0;
  double heu_lpg_budget_ms = 20.0;
  int64_t heu_lpg_node_interval = 200;
  int32_t heu_lpg_deterministic_restarts = 32;
  double heu_lpg_edge_threshold = 0.1;
  double heu_lpg_node_threshold = 0.5;
  double heu_lpg_lp_threshold = 0.1;
  double heu_lpg_seed_threshold = 0.3;

  // --- Heuristics: HiGHS sub-MIP ---
  bool heu_highs_submip_sec = true;

  // --- SPI separator ---
  int32_t spi_max_dp_size = 15;
  int32_t spi_max_seeds = 50;
  int32_t spi_max_grow_candidates = 40;

  // --- RC fixing ---
  RCFixingSettings rc_settings;

  // --- Labeling work limit ---
  int64_t labeling_max_queue_pops = 1'000'000;

  // --- Cut families ---
  bool enable_sec = true;
  bool enable_rci = true;
  bool enable_multistar = true;
  bool enable_comb = false;
  bool enable_rglm = false;
  bool enable_spi = false;
  int32_t max_cuts_sec = -1;
  int32_t max_cuts_rci = -1;
  int32_t max_cuts_multistar = -1;
  int32_t max_cuts_comb = -1;
  int32_t max_cuts_rglm = -1;
  int32_t max_cuts_spi = -1;

  // --- Other ---
  std::string branch_hyper = "off";
  double cutoff = std::numeric_limits<double>::infinity();
  int32_t preproc_fast_restarts = 12;
  int32_t hyper_sb_max_depth = 0;
  int32_t hyper_sb_iter_limit = 100;
  int32_t hyper_sb_min_reliable = 4;
  int32_t hyper_sb_max_candidates = 3;

  double user_time_limit = -1.0;  // -1 = no limit specified

  std::map<std::string, std::string> accepted_highs_options;
  bool threads_option_explicit = false;

  auto parse_bool = [](const std::string& v) {
    return v == "true" || v == "1";
  };

  for (const auto& [key, value] : options) {
    if (key == "presolve") {
      if (value != "off" && value != "false" && value != "0") {
        logger_.log(
            "Warning: option presolve={} requested, but cptp-bac currently "
            "forces presolve=off",
            value);
      }
      continue;
    }
    // --- Bounds & Preprocessing ---
    if (key == "labeling_elim_bounds") {
      labeling_elim_bounds = parse_bool(value);
      continue;
    }
    if (key == "all_pairs_bounds") {
      all_pairs_bounds = parse_bool(value);
      continue;
    }
    if (key == "edge_elimination") {
      edge_elimination = parse_bool(value);
      continue;
    }
    if (key == "edge_elimination_nodes") {
      edge_elimination_nodes = parse_bool(value);
      continue;
    }
    if (key == "bounds_propagation") {
      bounds_propagation = parse_bool(value);
      continue;
    }
    // --- Separation ---
    if (key == "separation_interval") {
      separation_interval = std::stoi(value);
      continue;
    }
    if (key == "separation_tol") {
      separation_tol = std::stod(value);
      continue;
    }
    if (key == "max_cuts_per_round") {
      max_cuts_per_round = std::stoi(value);
      continue;
    }
    // --- Heuristics: warm-start ---
    if (key == "heu_ws") {
      heu_ws = parse_bool(value);
      continue;
    }
    if (key == "heu_ws_ls_max_iter") {
      heu_ws_ls_max_iter = std::stoi(value);
      continue;
    }
    // --- Heuristics: LP-guided ---
    if (key == "heu_lpg") {
      heu_lpg = parse_bool(value);
      continue;
    }
    if (key == "heu_lpg_strategy") {
      heu_lpg_strategy = std::stoi(value);
      continue;
    }
    if (key == "heu_lpg_budget_ms") {
      heu_lpg_budget_ms = std::stod(value);
      continue;
    }
    if (key == "heu_lpg_node_interval") {
      heu_lpg_node_interval = std::stoll(value);
      continue;
    }
    if (key == "heu_lpg_deterministic_restarts") {
      heu_lpg_deterministic_restarts = std::stoi(value);
      continue;
    }
    if (key == "heu_lpg_edge_threshold") {
      heu_lpg_edge_threshold = std::stod(value);
      continue;
    }
    if (key == "heu_lpg_node_threshold") {
      heu_lpg_node_threshold = std::stod(value);
      continue;
    }
    if (key == "heu_lpg_lp_threshold") {
      heu_lpg_lp_threshold = std::stod(value);
      continue;
    }
    if (key == "heu_lpg_seed_threshold") {
      heu_lpg_seed_threshold = std::stod(value);
      continue;
    }
    // --- Heuristics: HiGHS sub-MIP ---
    if (key == "heu_highs_submip_sec") {
      heu_highs_submip_sec = parse_bool(value);
      continue;
    }
    // --- SPI separator ---
    if (key == "spi_max_dp_size") {
      spi_max_dp_size = std::stoi(value);
      continue;
    }
    if (key == "spi_max_seeds") {
      spi_max_seeds = std::stoi(value);
      continue;
    }
    if (key == "spi_max_grow_candidates") {
      spi_max_grow_candidates = std::stoi(value);
      continue;
    }
    // --- RC fixing ---
    if (key == "rc_fixing") {
      if (value == "root_only")
        rc_settings.strategy = RCFixingStrategy::root_only;
      else if (value == "on_ub_improvement")
        rc_settings.strategy = RCFixingStrategy::on_ub_improvement;
      else if (value == "periodic")
        rc_settings.strategy = RCFixingStrategy::periodic;
      else if (value == "adaptive")
        rc_settings.strategy = RCFixingStrategy::adaptive;
      else
        rc_settings.strategy = RCFixingStrategy::off;
      continue;
    }
    if (key == "rc_fixing_interval") {
      rc_settings.periodic_interval = std::stoi(value);
      continue;
    }
    if (key == "rc_fixing_to_one") {
      rc_settings.fix_to_one = parse_bool(value);
      continue;
    }
    if (key == "labeling_max_queue_pops") {
      labeling_max_queue_pops = std::stoll(value);
      continue;
    }
    // --- Cut families ---
    if (key == "enable_sec") {
      enable_sec = parse_bool(value);
      continue;
    }
    if (key == "enable_rci") {
      enable_rci = parse_bool(value);
      continue;
    }
    if (key == "enable_multistar") {
      enable_multistar = parse_bool(value);
      continue;
    }
    if (key == "enable_comb") {
      enable_comb = parse_bool(value);
      continue;
    }
    if (key == "enable_rglm") {
      enable_rglm = parse_bool(value);
      continue;
    }
    if (key == "enable_spi") {
      enable_spi = parse_bool(value);
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
    // --- Other ---
    if (key == "cutoff") {
      cutoff = std::stod(value);
      continue;
    }
    if (key == "branch_hyper") {
      branch_hyper = value;
      continue;
    }
    if (key == "preproc_fast_restarts") {
      preproc_fast_restarts = std::stoi(value);
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
    if (key == "output_flag") {
      output_flag = parse_bool(value);
      continue;
    }
    if (key == "random_seed") {
      random_seed = std::stoi(value);
      // Forward to HiGHS but don't duplicate in accepted_highs_options
      // (we log it ourselves in nondefault_settings).
      highs.setOptionValue(key, value);
      continue;
    }
    if (key == "time_limit") {
      user_time_limit = std::stod(value);
      // Don't set on HiGHS yet; we adjust for preprocessing time later.
      accepted_highs_options[key] = value;
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

  // Cascading disables: labeling_elim_bounds=false disables bounds-dependent
  // features
  if (!labeling_elim_bounds) {
    edge_elimination = false;
    edge_elimination_nodes = false;
    bounds_propagation = false;
    all_pairs_bounds = false;
    logger_.log(
        "labeling_elim_bounds=false: cascading off edge_elimination, "
        "edge_elimination_nodes, bounds_propagation, all_pairs_bounds");
  }

  preproc_fast_restarts = std::max<int32_t>(1, preproc_fast_restarts);
  heu_lpg_node_interval = std::max<int64_t>(1, heu_lpg_node_interval);
  heu_lpg_deterministic_restarts =
      std::max<int32_t>(1, heu_lpg_deterministic_restarts);
  heu_ws_ls_max_iter = std::max<int32_t>(1, heu_ws_ls_max_iter);
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
  const char* simd_str =
#ifdef __AVX512F__
      ", AVX-512";
#elif defined(__AVX2__)
      ", AVX2";
#elif defined(__AVX__)
      ", AVX";
#elif defined(__SSE4_2__)
      ", SSE4.2";
#else
      "";
#endif
  std::string highs_githash = highs.githash();
  if (highs_githash.empty()) highs_githash = "unknown";
  logger_.log("Running HiGHS {} (git hash: {})", highs.version(),
              highs_githash);
  logger_.log(
      "Thread count: {} physical cores, {} logical processors, using up to {} "
      "threads{}",
      physical_cores, logical_threads, preproc_thread_limit, simd_str);

  std::vector<std::string> nondefault_settings;
  if (random_seed != 0) {
    nondefault_settings.push_back("random_seed=" + std::to_string(random_seed));
  }
  nondefault_settings.push_back("presolve=off");
  if (!labeling_elim_bounds)
    nondefault_settings.push_back("labeling_elim_bounds=off");
  if (labeling_max_queue_pops != 1'000'000)
    nondefault_settings.push_back("labeling_max_queue_pops=" +
                                  std::to_string(labeling_max_queue_pops));
  if (all_pairs_bounds) nondefault_settings.push_back("all_pairs_bounds=on");
  if (!edge_elimination) nondefault_settings.push_back("edge_elimination=off");
  if (!edge_elimination_nodes)
    nondefault_settings.push_back("edge_elimination_nodes=off");
  if (!bounds_propagation)
    nondefault_settings.push_back("bounds_propagation=off");
  if (separation_interval != 1) {
    nondefault_settings.push_back("separation_interval=" +
                                  std::to_string(separation_interval));
  }
  if (max_cuts_per_round != 10) {
    nondefault_settings.push_back("max_cuts_per_round=" +
                                  std::to_string(max_cuts_per_round));
  }
  if (!heu_ws) nondefault_settings.push_back("heu_ws=off");
  if (heu_ws_ls_max_iter != 1000) {
    nondefault_settings.push_back("heu_ws_ls_max_iter=" +
                                  std::to_string(heu_ws_ls_max_iter));
  }
  if (!heu_lpg) nondefault_settings.push_back("heu_lpg=off");
  if (heu_lpg_strategy != 0) {
    nondefault_settings.push_back("heu_lpg_strategy=" +
                                  std::to_string(heu_lpg_strategy));
  }
  if (heu_lpg_node_interval != 200) {
    nondefault_settings.push_back("heu_lpg_node_interval=" +
                                  std::to_string(heu_lpg_node_interval));
  }
  if (heu_lpg_deterministic_restarts != 32) {
    nondefault_settings.push_back(
        "heu_lpg_deterministic_restarts=" +
        std::to_string(heu_lpg_deterministic_restarts));
  }
  if (heu_lpg_edge_threshold != 0.1)
    nondefault_settings.push_back("heu_lpg_edge_threshold=" +
                                  std::to_string(heu_lpg_edge_threshold));
  if (heu_lpg_node_threshold != 0.5)
    nondefault_settings.push_back("heu_lpg_node_threshold=" +
                                  std::to_string(heu_lpg_node_threshold));
  if (heu_lpg_lp_threshold != 0.1)
    nondefault_settings.push_back("heu_lpg_lp_threshold=" +
                                  std::to_string(heu_lpg_lp_threshold));
  if (heu_lpg_seed_threshold != 0.3)
    nondefault_settings.push_back("heu_lpg_seed_threshold=" +
                                  std::to_string(heu_lpg_seed_threshold));
  if (!heu_highs_submip_sec)
    nondefault_settings.push_back("heu_highs_submip_sec=off");
  if (spi_max_dp_size != 15)
    nondefault_settings.push_back("spi_max_dp_size=" +
                                  std::to_string(spi_max_dp_size));
  if (spi_max_seeds != 50)
    nondefault_settings.push_back("spi_max_seeds=" +
                                  std::to_string(spi_max_seeds));
  if (spi_max_grow_candidates != 40)
    nondefault_settings.push_back("spi_max_grow_candidates=" +
                                  std::to_string(spi_max_grow_candidates));
  if (branch_hyper != "off")
    nondefault_settings.push_back("branch_hyper=" + branch_hyper);
  if (preproc_fast_restarts != 12) {
    nondefault_settings.push_back("preproc_fast_restarts=" +
                                  std::to_string(preproc_fast_restarts));
  }
  if (cutoff < std::numeric_limits<double>::infinity()) {
    nondefault_settings.push_back("cutoff=" + std::to_string(cutoff));
  }
  if (heu_lpg_budget_ms != 20.0) {
    nondefault_settings.push_back("heu_lpg_budget_ms=" +
                                  std::to_string(heu_lpg_budget_ms));
  }
  if (separation_tol != sep::kDefaultFracTol) {
    nondefault_settings.push_back("separation_tol=" +
                                  std::to_string(separation_tol));
  }
  if (rc_settings.strategy != RCFixingStrategy::adaptive) {
    auto rc_str = [](RCFixingStrategy s) -> const char* {
      switch (s) {
        case RCFixingStrategy::root_only:
          return "root_only";
        case RCFixingStrategy::on_ub_improvement:
          return "on_ub_improvement";
        case RCFixingStrategy::periodic:
          return "periodic";
        case RCFixingStrategy::adaptive:
          return "adaptive";
        default:
          return "off";
      }
    };
    nondefault_settings.push_back(std::string("rc_fixing=") +
                                  rc_str(rc_settings.strategy));
    if (rc_settings.periodic_interval != 100) {
      nondefault_settings.push_back(
          "rc_fixing_interval=" +
          std::to_string(rc_settings.periodic_interval));
    }
    if (rc_settings.fix_to_one) {
      nondefault_settings.push_back("rc_fixing_to_one=on");
    }
  }
  if (!enable_sec) nondefault_settings.push_back("enable_sec=off");
  if (!enable_rci) nondefault_settings.push_back("enable_rci=off");
  if (!enable_multistar) nondefault_settings.push_back("enable_multistar=off");
  if (enable_comb) nondefault_settings.push_back("enable_comb=on");
  if (enable_spi) nondefault_settings.push_back("enable_spi=on");
  if (enable_rglm) nondefault_settings.push_back("enable_rglm=on");
  if (max_cuts_sec != -1)
    nondefault_settings.push_back("max_cuts_sec=" +
                                  std::to_string(max_cuts_sec));
  if (max_cuts_rci != -1)
    nondefault_settings.push_back("max_cuts_rci=" +
                                  std::to_string(max_cuts_rci));
  if (max_cuts_multistar != -1)
    nondefault_settings.push_back("max_cuts_multistar=" +
                                  std::to_string(max_cuts_multistar));
  if (max_cuts_comb != -1)
    nondefault_settings.push_back("max_cuts_comb=" +
                                  std::to_string(max_cuts_comb));
  if (max_cuts_rglm != -1)
    nondefault_settings.push_back("max_cuts_rglm=" +
                                  std::to_string(max_cuts_rglm));
  if (max_cuts_spi != -1)
    nondefault_settings.push_back("max_cuts_spi=" +
                                  std::to_string(max_cuts_spi));
  if (branch_hyper != "off") {
    if (hyper_sb_max_depth != 0)
      nondefault_settings.push_back("branch_hyper_sb_max_depth=" +
                                    std::to_string(hyper_sb_max_depth));
    if (hyper_sb_iter_limit != 100)
      nondefault_settings.push_back("branch_hyper_sb_iter_limit=" +
                                    std::to_string(hyper_sb_iter_limit));
    if (hyper_sb_min_reliable != 4)
      nondefault_settings.push_back("branch_hyper_sb_min_reliable=" +
                                    std::to_string(hyper_sb_min_reliable));
    if (hyper_sb_max_candidates != 3)
      nondefault_settings.push_back("branch_hyper_sb_max_candidates=" +
                                    std::to_string(hyper_sb_max_candidates));
  }
  for (const auto& [k, v] : accepted_highs_options) {
    nondefault_settings.push_back(k + "=" + v);
  }

  const int32_t source = problem_.source();
  const int32_t target = problem_.target();
  const int32_t n = problem_.num_nodes();
  const int32_t m = problem_.num_edges();
  const std::string settings_summary = join_tokens(nondefault_settings);
  const std::string solve_label =
      problem_.name.empty() ? std::string("instance") : problem_.name;
  logger_.log("Solving {} {} with:", solve_label,
              problem_.is_tour() ? "tour problem" : "path problem");
  if (problem_.is_tour()) {
    logger_.log("  {} nodes, {} edges, tour through node {}", n, m, source);
  } else {
    logger_.log("  {} nodes, {} edges, path from node {} to node {}", n, m,
                source, target);
  }
  logger_.log("  Settings  : {}", settings_summary);

  // When cutoff is provided, set HiGHS objective_bound for node pruning.
  // This sets the initial upper_limit in the MIP solver, enabling pruning
  // of nodes whose LP bound exceeds the cutoff.
  // When heuristics are disabled, turn off HiGHS internal MIP heuristics.
  if (cutoff < std::numeric_limits<double>::infinity()) {
    highs.setOptionValue("objective_bound", cutoff);
  }
  // Note: HiGHS internal heuristics (mip_heuristic_effort) are controlled
  // separately via the HiGHS option passthrough, not by heu_ws.
  // Enforce disabled presolve even if user passed --presolve.
  highs.setOptionValue("presolve", "off");

  // Honor user output preference for HiGHS logging (single-source output).
  highs.setOptionValue("output_flag", output_flag);
  if (output_flag) {
#ifdef _WIN32
    constexpr const char* kNullLogPath = "NUL";
#else
    constexpr const char* kNullLogPath = "/dev/null";
#endif
    // Route HiGHS logging through callback so we can filter specific lines.
    highs.setOptionValue("log_file", kNullLogPath);
    highs.setOptionValue("log_to_console", false);
  }

  // Deterministic staged preprocessing:
  // 1) Stage-1 bounds (parallel forward/backward)
  // 2) Construction seeds (1/2-customer)
  // 3) Local search on top seeds (parallel)
  // 3.5) Edge elimination from construction UB
  // 4) Reduced local search + deeper ng/DSSR (parallel)
  // 5) Start HiGHS
  // correction = profit(source) when s == t (tour: depot profit
  // double-subtracted)
  double correction = problem_.is_tour() ? problem_.profit(source) : 0.0;

  std::vector<double> fwd_bounds, bwd_bounds;
  std::vector<double> all_pairs;
  heuristic::HeuristicResult warm_start{
      {}, std::numeric_limits<double>::infinity()};
  double warm_start_ub = std::numeric_limits<double>::infinity();
  int32_t stage3_elim_edges = 0;
  double stage3_elim_ratio = 0.0;
  int32_t stage5_elim_edges = 0;
  double stage5_elim_ratio = 0.0;

  // Use cutoff as UB when provided (known optimal for prove-only mode)
  if (cutoff < std::numeric_limits<double>::infinity()) {
    warm_start_ub = cutoff;
  }

  {
    Timer stage1_timer;
    if (labeling_elim_bounds) {
      logger_.log("Lower bounds calculation:");
      logger_.log("  capacity-aware labeling (source={}, target={})", source,
                  target);
      if (source != target) {
        std::jthread t_fwd([&] {
          fwd_bounds = preprocess::labeling_from(problem_, source,
                                                 labeling_max_queue_pops);
        });
        bwd_bounds =
            preprocess::labeling_from(problem_, target, labeling_max_queue_pops);
      } else {
        fwd_bounds = preprocess::labeling_from(problem_, source,
                                               labeling_max_queue_pops);
        bwd_bounds = fwd_bounds;
      }
      if (fwd_bounds.empty() || bwd_bounds.empty()) {
        logger_.log(
            "  labeling budget exhausted ({} pops) — disabling "
            "bounds-dependent features",
            labeling_max_queue_pops);
        fwd_bounds.clear();
        bwd_bounds.clear();
        labeling_elim_bounds = false;
        edge_elimination = false;
        edge_elimination_nodes = false;
        bounds_propagation = false;
        all_pairs_bounds = false;
        rc_settings.strategy = RCFixingStrategy::off;
      } else {
        logger_.log(
            "  fwd_reachable={}/{}  bwd_reachable={}/{}  fwd[target]={}  "
            "time={:.3f}s",
            count_finite_bounds(fwd_bounds), n, count_finite_bounds(bwd_bounds),
            n,
            (target >= 0 && target < n &&
             target < static_cast<int32_t>(fwd_bounds.size()))
                ? fwd_bounds[static_cast<size_t>(target)]
                : std::numeric_limits<double>::infinity(),
            stage1_timer.elapsed_seconds());
      }
    } else {
      // No log when bounds skipped (matches PR #49 style)
    }
    heuristic::ConstructionPool construction_pool;
    heuristic::HeuristicResult construction_start{
        {}, std::numeric_limits<double>::infinity()};
    heuristic::HeuristicResult ls_stage4{
        {}, std::numeric_limits<double>::infinity()};
    double construction_ub = std::numeric_limits<double>::infinity();
    std::vector<bool> stage3_edge_active;
    logger_.log("Construction heuristic:");
    Timer stage2_timer;

    if (heu_ws) {
      const int32_t candidate_cap =
          std::max<int32_t>(preproc_fast_restarts, preproc_thread_limit);
      construction_pool =
          heuristic::build_construction_pool(problem_, candidate_cap);
      construction_start =
          heuristic::best_construction_solution(problem_, construction_pool);
      if (has_feasible_heuristic(construction_start)) {
        construction_ub = construction_start.objective;
        warm_start = construction_start;
        logger_.log("  ub={}  candidates={}  time={:.3f}s", construction_ub,
                    static_cast<int32_t>(construction_pool.candidates.size()),
                    stage2_timer.elapsed_seconds());
      } else {
        logger_.log("  no feasible incumbent  time={:.3f}s",
                    stage2_timer.elapsed_seconds());
      }
    } else {
      logger_.log("  skipped (heu_ws=false)  time={:.3f}s",
                  stage2_timer.elapsed_seconds());
    }

    if (cutoff >= std::numeric_limits<double>::infinity()) {
      warm_start_ub = has_feasible_heuristic(warm_start)
                          ? warm_start.objective
                          : std::numeric_limits<double>::infinity();
    } else if (has_feasible_heuristic(warm_start)) {
      warm_start_ub = std::min(warm_start_ub, warm_start.objective);
    }

    logger_.log("Preprocess:");
    Timer stage3_timer;
    if (edge_elimination && std::isfinite(construction_ub) &&
        !fwd_bounds.empty() && !bwd_bounds.empty()) {
      auto eliminated = preprocess::edge_elimination(
          problem_, fwd_bounds, bwd_bounds, construction_ub, correction);
      stage3_elim_edges = static_cast<int32_t>(
          std::count(eliminated.begin(), eliminated.end(), true));
      stage3_elim_ratio = (m > 0) ? static_cast<double>(stage3_elim_edges) /
                                        static_cast<double>(m)
                                  : 0.0;
      stage3_edge_active.assign(static_cast<size_t>(m), true);
      for (int32_t e = 0; e < m; ++e) {
        if (eliminated[static_cast<size_t>(e)]) {
          stage3_edge_active[static_cast<size_t>(e)] = false;
        }
      }
      logger_.log("  eliminated={}/{} ({:.1f}%)  time={:.3f}s",
                  stage3_elim_edges, m, 100.0 * stage3_elim_ratio,
                  stage3_timer.elapsed_seconds());
    } else {
      logger_.log(
          "  skipped (edge_elimination={}, finite_ub={}, have_bounds={})  "
          "time={:.3f}s",
          edge_elimination ? "true" : "false",
          std::isfinite(construction_ub) ? "true" : "false",
          (!fwd_bounds.empty() && !bwd_bounds.empty()) ? "true" : "false",
          stage3_timer.elapsed_seconds());
    }

    if (all_pairs_bounds) {
      logger_.log("Lower bounds all pairs calculation:");
      Timer all_pairs_timer;
      {
        constexpr double inf = std::numeric_limits<double>::infinity();
        all_pairs.assign(static_cast<size_t>(n) * n, inf);
        parallel::parallel_for(
            0, n,
            [&](int s) {
              auto row = preprocess::forward_labeling(problem_, s);
              std::copy(row.begin(), row.end(),
                        all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
            },
            preproc_thread_limit);
      }
      if (!all_pairs.empty()) {
        const int64_t all_pairs_total =
            static_cast<int64_t>(n) * static_cast<int64_t>(n);
        logger_.log("  all_pairs_reachable={}/{}  time={:.3f}s",
                    count_finite_bounds(all_pairs), all_pairs_total,
                    all_pairs_timer.elapsed_seconds());
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
        logger_.log("  no bounds produced");
      }
    }

    logger_.log("Local search:");
    Timer stage4_timer;
    if (heu_ws && !construction_pool.candidates.empty()) {
      const int64_t scaled_starts =
          static_cast<int64_t>(preproc_thread_limit) *
          static_cast<int64_t>(kWarmStartLsStartsPerThread);
      const int32_t starts_from_threads =
          static_cast<int32_t>(std::clamp<int64_t>(
              scaled_starts, 1, std::numeric_limits<int32_t>::max()));
      const int32_t requested_starts =
          std::max<int32_t>(preproc_fast_restarts, starts_from_threads);
      logger_.log("  starts={}  max_iter_per_start={}", requested_starts,
                  heu_ws_ls_max_iter);
      logger_.log("  starts iter_accum           ub impr     time");
      heuristic::WarmStartProgressOptions progress_opts;
      progress_opts.report_every_starts = 32;
      progress_opts.report_on_ub_improvement = true;
      auto progress_cb = [&](const heuristic::WarmStartProgressSnapshot& snap) {
        logger_.log("  {:>6} {:>10} {:>12.6g} {:>4} {:>8.3f}s",
                    snap.starts_done, snap.ls_iterations_total,
                    snap.best_objective, snap.ub_improvements,
                    stage4_timer.elapsed_seconds());
      };
      ls_stage4 = heuristic::run_local_search_from_pool(
          problem_, construction_pool, requested_starts, heu_ws_ls_max_iter,
          preproc_thread_limit, stage3_edge_active, &progress_opts, progress_cb,
          true);
    } else {
      logger_.log("  skipped (heu_ws=false or candidates=0)  time={:.3f}s",
                  stage4_timer.elapsed_seconds());
    }

    if (has_feasible_heuristic(ls_stage4) &&
        (!has_feasible_heuristic(warm_start) ||
         ls_stage4.objective + 1e-9 < warm_start.objective)) {
      warm_start = std::move(ls_stage4);
    }

    if (cutoff >= std::numeric_limits<double>::infinity()) {
      warm_start_ub = has_feasible_heuristic(warm_start)
                          ? warm_start.objective
                          : std::numeric_limits<double>::infinity();
    } else if (has_feasible_heuristic(warm_start)) {
      warm_start_ub = std::min(warm_start_ub, warm_start.objective);
    }
    logger_.log("Preprocess restart:");
    Timer stage5_timer;
    if (edge_elimination && std::isfinite(warm_start_ub) &&
        warm_start_ub < construction_ub - 1e-9 && !fwd_bounds.empty() &&
        !bwd_bounds.empty()) {
      auto eliminated = preprocess::edge_elimination(
          problem_, fwd_bounds, bwd_bounds, warm_start_ub, correction);
      stage5_elim_edges = static_cast<int32_t>(
          std::count(eliminated.begin(), eliminated.end(), true));
      stage5_elim_ratio = (m > 0) ? static_cast<double>(stage5_elim_edges) /
                                        static_cast<double>(m)
                                  : 0.0;
      logger_.log("  eliminated={}/{} ({:.1f}%)  time={:.3f}s",
                  stage5_elim_edges, m, 100.0 * stage5_elim_ratio,
                  stage5_timer.elapsed_seconds());
    } else if (edge_elimination && std::isfinite(warm_start_ub)) {
      logger_.log("  skipped (ub not improved)  time=0.000s");
    } else {
      logger_.log(
          "  skipped (edge_elimination={}, finite_ub={}, have_bounds={})  "
          "time={:.3f}s",
          edge_elimination ? "true" : "false",
          std::isfinite(warm_start_ub) ? "true" : "false",
          (!fwd_bounds.empty() && !bwd_bounds.empty()) ? "true" : "false",
          stage5_timer.elapsed_seconds());
    }
  }
  std::vector<int8_t> fixed_y(static_cast<size_t>(n), static_cast<int8_t>(-1));
  HiGHSBridge bridge(problem_, highs, logger_, separation_tol);
  bridge.set_separation_interval(separation_interval);
  bridge.set_max_cuts_per_round(max_cuts_per_round);
  bridge.set_submip_separation(heu_highs_submip_sec);
  bridge.set_upper_bound(warm_start_ub);
  bridge.set_rc_fixing(rc_settings);
  bridge.set_labeling_max_queue_pops(labeling_max_queue_pops);
  bridge.set_edge_elimination(edge_elimination);
  bridge.set_edge_elimination_nodes(edge_elimination_nodes);
  bridge.set_heuristic_node_interval(heu_lpg_node_interval);
  bridge.set_heuristic_deterministic_restarts(heu_lpg_deterministic_restarts);
  bridge.set_fixed_y(fixed_y);
  if (max_cuts_sec >= 0) bridge.set_separator_max_cuts("SEC", max_cuts_sec);
  if (max_cuts_rci >= 0) bridge.set_separator_max_cuts("RCI", max_cuts_rci);
  if (max_cuts_multistar >= 0)
    bridge.set_separator_max_cuts("Multistar", max_cuts_multistar);
  if (max_cuts_comb >= 0) bridge.set_separator_max_cuts("Comb", max_cuts_comb);
  if (max_cuts_rglm >= 0) bridge.set_separator_max_cuts("RGLM", max_cuts_rglm);
  if (max_cuts_spi >= 0) bridge.set_separator_max_cuts("SPI", max_cuts_spi);

  bridge.set_labeling_bounds(std::move(fwd_bounds), std::move(bwd_bounds),
                             correction);
  if (all_pairs_bounds) {
    bridge.set_all_pairs_bounds(std::move(all_pairs));
  }

  // Add separators (same default families for tour and s-t path).
  if (enable_sec) bridge.add_separator(std::make_unique<sep::SECSeparator>());
  if (enable_rci) bridge.add_separator(std::make_unique<sep::RCISeparator>());
  if (enable_multistar)
    bridge.add_separator(std::make_unique<sep::MultistarSeparator>());
  if (enable_rglm) bridge.add_separator(std::make_unique<sep::RGLMSeparator>());
  if (enable_comb) bridge.add_separator(std::make_unique<sep::CombSeparator>());
  if (all_pairs_bounds && enable_spi) {
    auto spi = std::make_unique<sep::SPISeparator>();
    spi->set_max_dp_size(spi_max_dp_size);
    spi->set_max_seeds(spi_max_seeds);
    spi->set_max_grow_candidates(spi_max_grow_candidates);
    bridge.add_separator(std::move(spi));
  }

  bridge.build_formulation();

  configure_hyperplane_branching(
      problem_, bridge, logger_, branch_hyper, hyper_sb_max_depth,
      hyper_sb_iter_limit, hyper_sb_min_reliable, hyper_sb_max_candidates);

  bridge.install_separators();
  bridge.set_bounds_propagation(bounds_propagation);
  bridge.install_propagator();
  bridge.set_heuristic_callback(heu_lpg);
  bridge.set_heuristic_budget_ms(heu_lpg_budget_ms);
  bridge.set_heuristic_strategy(heu_lpg_strategy);
  bridge.set_heu_lpg_edge_threshold(heu_lpg_edge_threshold);
  bridge.set_heu_lpg_node_threshold(heu_lpg_node_threshold);
  bridge.set_heu_lpg_lp_threshold(heu_lpg_lp_threshold);
  bridge.set_heu_lpg_seed_threshold(heu_lpg_seed_threshold);
  bridge.install_heuristic_callback();

  if (heu_ws && has_feasible_heuristic(warm_start)) {
    // Pass heuristic solution to HiGHS (skip if heuristics disabled)
    HighsSolution start;
    start.value_valid = true;
    start.col_value = std::move(warm_start.col_values);
    HighsInt num_cols = highs.getNumCol();
    while (static_cast<HighsInt>(start.col_value.size()) < num_cols)
      start.col_value.push_back(0.0);
    highs.setSolution(start);
  }

  const double pre_highs_seconds = timer.elapsed_seconds();
  // Adjust time_limit to account for preprocessing time already spent.
  if (user_time_limit > 0.0) {
    double remaining = user_time_limit - pre_highs_seconds;
    if (remaining < 0.1) remaining = 0.1;  // give HiGHS at least 0.1s
    highs.setOptionValue("time_limit", remaining);
  }
  logger_.log("HiGHS start offset (preprocessing): {:.3f}s", pre_highs_seconds);
  highs.run();
  const double total_elapsed_seconds = timer.elapsed_seconds();
  logger_.log("Total elapsed (preprocessing + HiGHS): {:.3f}s",
              total_elapsed_seconds);

  // Clear static callbacks/state to avoid leaking between solves
  HighsUserSeparator::clearCallback();

  auto result = bridge.extract_result();
  result.time_seconds = total_elapsed_seconds;

  return result;
}

}  // namespace cptp
