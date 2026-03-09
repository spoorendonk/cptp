#include "model/highs_bridge.h"

#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <mutex>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "heuristic/primal_heuristic.h"
#include "lp_data/HighsCallback.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsUserPropagator.h"
#include "mip/HighsUserSeparator.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/reachability.h"
#include "sep/sec_separator.h"
#include "sep/separation_context.h"
#include "sep/separation_oracle.h"

namespace cptp {

namespace {

struct PropagatorAdjEntry {
  int32_t edge;
  int32_t neighbor;
};

bool should_run_rc_fixing(const RCFixingSettings& settings, bool ub_improved,
                          int64_t ub_improvements, int64_t propagator_calls,
                          bool adaptive_disabled) {
  switch (settings.strategy) {
    case RCFixingStrategy::root_only:
      // Run once when UB first becomes available.
      return ub_improved && ub_improvements == 1;
    case RCFixingStrategy::on_ub_improvement:
      return ub_improved;
    case RCFixingStrategy::periodic:
      return (propagator_calls % settings.periodic_interval) == 0;
    case RCFixingStrategy::adaptive:
      return ub_improved && !adaptive_disabled;
    case RCFixingStrategy::off:
    default:
      return false;
  }
}

}  // namespace

HiGHSBridge::HiGHSBridge(const Problem& prob, Highs& highs, Logger& logger,
                         double frac_tol)
    : prob_(prob),
      highs_(highs),
      logger_(logger),
      num_edges_(prob.num_edges()),
      num_nodes_(prob.num_nodes()),
      frac_tol_(frac_tol),
      int_tol_(1e-6),
      cached_x_lp_(std::make_shared<std::vector<double>>()),
      cached_y_lp_(std::make_shared<std::vector<double>>()),
      lp_cache_mutex_(std::make_shared<std::mutex>()),
      heuristic_calls_(std::make_shared<int64_t>(0)),
      heuristic_solutions_(std::make_shared<int64_t>(0)),
      heuristic_time_seconds_(std::make_shared<double>(0.0)) {
  // Integral feasibility: tight but not zero to avoid rejecting
  // valid solutions due to floating-point noise in LP values.
  // Fractional separation uses a looser tolerance (default 1e-2).
}

HiGHSBridge::~HiGHSBridge() {
  HighsUserSeparator::clearCallback();
  HighsUserPropagator::clearCallback();
}

void HiGHSBridge::add_separator(std::unique_ptr<sep::Separator> sep) {
  if (sep) separator_stats_.try_emplace(sep->name(), SeparatorStats{});
  separators_.push_back(std::move(sep));
}

void HiGHSBridge::build_formulation() {
  const auto& graph = prob_.graph();
  const int32_t n = num_nodes_;
  const int32_t m = num_edges_;
  const int total_vars = m + n;

  std::vector<double> col_cost(total_vars, 0.0);
  std::vector<double> col_lower(total_vars, 0.0);
  std::vector<double> col_upper(total_vars, 1.0);

  // Tour mode: depot-incident edges allow x_e ∈ {0,1,2} to permit 2-node tours
  // (Jepsen et al. 2014, standard CPTP/CVRP formulation)
  if (prob_.is_tour()) {
    for (auto e : graph.incident_edges(prob_.source())) {
      col_upper[e] = 2.0;
    }
  }

  // Fix source and target y variables to 1 via bounds
  col_lower[m + prob_.source()] = 1.0;
  col_lower[m + prob_.target()] = 1.0;

  // Fix unreachable nodes: demand-reachability preprocessing
  if (prob_.capacity() > 0 && prob_.capacity() < 1e17) {
    auto reachable = preprocess::demand_reachability(prob_);
    for (int32_t i = 0; i < n; ++i) {
      if (!reachable[i] && i != prob_.source() && i != prob_.target()) {
        col_upper[m + i] = 0.0;
        for (auto e : graph.incident_edges(i)) col_upper[e] = 0.0;
      }
    }
  }

  // Optional fixed-y assignments supplied by the caller.
  if (!fixed_y_.empty() && fixed_y_.size() == static_cast<size_t>(n)) {
    for (int32_t i = 0; i < n; ++i) {
      if (i == prob_.source() || i == prob_.target()) continue;
      const int8_t fix = fixed_y_[static_cast<size_t>(i)];
      if (fix == 0) {
        col_upper[m + i] = 0.0;
        for (auto e : graph.incident_edges(i)) {
          col_upper[e] = 0.0;
        }
      } else if (fix == 1) {
        col_lower[m + i] = 1.0;
      }
    }
  }

  // Edge elimination: capacity-aware labeling bounds (staged preprocessing).
  if (edge_elimination_enabled_ &&
      upper_bound_ < std::numeric_limits<double>::infinity() &&
      !fwd_bounds_.empty() && !bwd_bounds_.empty()) {
    auto eliminated = preprocess::edge_elimination(
        prob_, fwd_bounds_, bwd_bounds_, upper_bound_, correction_);
    int32_t edge_elim_count = 0;
    for (auto e : graph.edges()) {
      if (eliminated[e] && col_upper[e] > 0.0) {
        col_upper[e] = 0.0;
        edge_elim_count++;
      }
    }
    // Also eliminate nodes where all incident edges are eliminated
    int32_t node_elim_count = 0;
    if (edge_elimination_nodes_) {
      for (int32_t i = 0; i < n; ++i) {
        if (i == prob_.source() || i == prob_.target() ||
            col_upper[m + i] == 0.0)
          continue;
        bool all_eliminated = true;
        for (auto e : graph.incident_edges(i)) {
          if (col_upper[e] > 0.0) {
            all_eliminated = false;
            break;
          }
        }
        if (all_eliminated) {
          col_upper[m + i] = 0.0;
          node_elim_count++;
        }
      }
    }
    if (edge_elim_count > 0 || node_elim_count > 0) {
      logger_.log("Edge elimination: {}/{} edges, {}/{} nodes fixed to 0",
                  edge_elim_count, m, node_elim_count, n - 1);
    }
  }

  // Objective: max sum(profit_i * y_i) - sum(cost_e * x_e)
  // HiGHS minimizes: min sum(cost_e * x_e) - sum(profit_i * y_i)
  for (auto e : graph.edges()) {
    col_cost[e] = prob_.edge_cost(e);
  }
  for (int32_t i = 0; i < n; ++i) {
    col_cost[m + i] = -prob_.profit(i);
  }

  // Check if all objective coefficients are integer-valued.
  obj_is_integer_ = std::all_of(col_cost.begin(), col_cost.end(), [](double c) {
    return std::abs(c - std::round(c)) < 1e-9;
  });

  // Add columns with costs
  for (int i = 0; i < total_vars; ++i) {
    highs_.addCol(col_cost[i], col_lower[i], col_upper[i], 0, nullptr, nullptr);
  }

  // Set integer types
  for (int i = 0; i < total_vars; ++i) {
    highs_.changeColIntegrality(i, HighsVarType::kInteger);
  }

  highs_.changeObjectiveSense(ObjSense::kMinimize);

  // 1. Degree constraints
  // Tour (source == target): sum x_e = 2*y_i for all nodes
  // Path (source != target): sum x_e = y_i for source/target (degree 1),
  //                          sum x_e = 2*y_i for intermediates (degree 2)
  for (int32_t i = 0; i < n; ++i) {
    std::vector<int> indices;
    std::vector<double> values;
    for (auto e : graph.incident_edges(i)) {
      indices.push_back(static_cast<int>(e));
      values.push_back(1.0);
    }
    bool is_terminal =
        !prob_.is_tour() && (i == prob_.source() || i == prob_.target());
    indices.push_back(m + i);
    values.push_back(is_terminal ? -1.0 : -2.0);
    highs_.addRow(0.0, 0.0, static_cast<int>(indices.size()), indices.data(),
                  values.data());
  }

  // 2. Capacity constraint: sum(demand_i * y_i) <= Q
  if (prob_.capacity() > 0 && prob_.capacity() < 1e17) {
    std::vector<int> indices;
    std::vector<double> values;
    for (int32_t i = 0; i < n; ++i) {
      double d = prob_.demand(i);
      if (d > 0) {
        indices.push_back(m + i);
        values.push_back(d);
      }
    }
    if (!indices.empty()) {
      highs_.addRow(-kHighsInf, prob_.capacity(),
                    static_cast<int>(indices.size()), indices.data(),
                    values.data());
    }
  }
}

void HiGHSBridge::install_separators() {
  const bool sec_enabled_for_callbacks =
      std::any_of(separators_.begin(), separators_.end(),
                  [](const auto& s) { return s && s->name() == "SEC"; });
  // Capture what the callback needs by pointer/reference.
  // The bridge must outlive the highs.run() call.
  HighsUserSeparator::setCallback([this, sec_enabled_for_callbacks](
                                      const HighsLpRelaxation& lpRelaxation,
                                      HighsCutPool& cutpool,
                                      const HighsMipSolver& mipsolver) {
    // Sub-MIPs have a different (presolved) column space.
    // When enabled, expand LP solution via undoPrimal and translate
    // cuts back to reduced space. Only run at the sub-MIP root node
    // to tighten the LP without per-node overhead.
    if (mipsolver.submip) {
      if (!submip_separation_) return;
      if (mipsolver.mipdata_->num_nodes > 0) return;  // root only
    }

    // Amortized separation: skip rounds to reduce overhead
    // (not applied for sub-MIPs — they are short-lived)
    if (!mipsolver.submip) {
      int64_t call = separation_calls_++;
      if (separation_interval_ > 1 && (call % separation_interval_) != 0)
        return;
    }

    const auto& sol = lpRelaxation.getSolution();
    const int32_t m = num_edges_;
    const int32_t n = num_nodes_;

    std::vector<double> x_vals(m);
    std::vector<double> y_vals(n);
    std::vector<int> orig_to_reduced;  // empty for main MIP

    if (mipsolver.submip) {
      // Sub-MIP: LP solution is in reduced (presolved) column space.
      // Expand to original space so separators see the full model.
      auto& postSolveStack = const_cast<presolve::HighsPostsolveStack&>(
          mipsolver.mipdata_->postSolveStack);

      HighsSolution temp_sol;
      temp_sol.col_value = sol.col_value;
      temp_sol.value_valid = true;
      postSolveStack.undoPrimal(*mipsolver.options_mip_, temp_sol);

      if (static_cast<int32_t>(temp_sol.col_value.size()) < m + n) return;

      std::copy_n(temp_sol.col_value.begin(), m, x_vals.begin());
      std::copy_n(temp_sol.col_value.begin() + m, n, y_vals.begin());

      // Build inverse column mapping (original → reduced) for cut translation
      const auto& pss = mipsolver.mipdata_->postSolveStack;
      int orig_ncol = pss.getOrigNumCol();
      int reduced_ncol = static_cast<int>(sol.col_value.size());
      orig_to_reduced.assign(orig_ncol, -1);
      for (int r = 0; r < reduced_ncol; r++) {
        orig_to_reduced[pss.getOrigColIndex(r)] = r;
      }
    } else {
      std::copy_n(sol.col_value.begin(), m, x_vals.begin());
      std::copy_n(sol.col_value.begin() + m, n, y_vals.begin());

      // Cache LP values for heuristic callback
      {
        std::lock_guard lock(*lp_cache_mutex_);
        *cached_x_lp_ = x_vals;
        *cached_y_lp_ = y_vals;
      }
    }

    // Build support graph once for all separators
    const auto& graph = prob_.graph();
    const double graph_tol = 1e-6;  // for edge inclusion in support graph
    digraph_builder builder(n);
    for (auto e : graph.edges()) {
      double xval = x_vals[e];
      if (xval > graph_tol) {
        int32_t u = graph.edge_source(e);
        int32_t v = graph.edge_target(e);
        builder.add_arc(u, v, xval);
        builder.add_arc(v, u, xval);
      }
    }
    auto [support_graph, capacity] = builder.build();

    // Build Gomory-Hu tree: n-1 max-flows instead of 3*n
    gomory_hu_tree flow_tree(support_graph, capacity, prob_.source());

    sep::SeparationContext ctx{
        .problem = prob_,
        .x_values = std::span(x_vals.data(), static_cast<size_t>(m)),
        .y_values = std::span(y_vals.data(), static_cast<size_t>(n)),
        .x_offset = x_offset(),
        .y_offset = y_offset(),
        .tol = frac_tol_,
        .flow_tree = &flow_tree,
        .upper_bound = mipsolver.mipdata_->upper_limit,
        .all_pairs = all_pairs_.empty() ? std::span<const double>{}
                                        : std::span<const double>(all_pairs_),
    };

    const bool is_submip = !orig_to_reduced.empty();
    if (is_submip) {
      if (!sec_enabled_for_callbacks) return;
      // Sub-MIP: only run SEC (lazy constraint) to keep overhead low.
      // Other separators (RCI, Multistar, etc.) add tightening cuts
      // but aren't needed for feasibility in short-lived sub-MIPs.
      sep::SECSeparator sec;
      auto cuts = sec.separate(ctx);
      const auto it_cap = per_separator_max_cuts_.find(sec.name());
      const int32_t max_cuts = (it_cap != per_separator_max_cuts_.end())
                                   ? it_cap->second
                                   : max_cuts_per_sep_;

      std::stable_sort(cuts.begin(), cuts.end(),
                       [](const sep::Cut& a, const sep::Cut& b) {
                         return a.violation > b.violation;
                       });
      int32_t added = 0;
      if (max_cuts == 0) return;
      for (int32_t j = 0; j < static_cast<int32_t>(cuts.size()); ++j) {
        if (max_cuts >= 0 && added >= max_cuts) break;
        auto& cut = cuts[j];
        // Translate cut indices from original to reduced space
        std::vector<HighsInt> red_idx;
        std::vector<double> red_val;
        double rhs = cut.rhs;

        for (size_t k = 0; k < cut.indices.size(); k++) {
          int orig_col = cut.indices[k];
          int red_col = orig_to_reduced[orig_col];
          if (red_col >= 0) {
            red_idx.push_back(red_col);
            red_val.push_back(cut.values[k]);
          } else {
            // Column eliminated by presolve — adjust RHS
            double fixed_val =
                (orig_col < m) ? x_vals[orig_col] : y_vals[orig_col - m];
            rhs -= cut.values[k] * fixed_val;
          }
        }

        if (!red_idx.empty()) {
          cutpool.addCut(mipsolver, red_idx.data(), red_val.data(),
                         static_cast<HighsInt>(red_idx.size()), rhs);
          total_cuts_++;
          added++;
        }
      }
    } else {
      // Main MIP: run all separators in parallel, timing each one
      using clock = std::chrono::steady_clock;
      std::vector<std::vector<sep::Cut>> results(separators_.size());
      std::vector<double> elapsed(separators_.size(), 0.0);
      tbb::task_group tg;
      for (size_t i = 0; i < separators_.size(); ++i) {
        tg.run([&, i] {
          auto t0 = clock::now();
          results[i] = separators_[i]->separate(ctx);
          elapsed[i] = std::chrono::duration<double>(clock::now() - t0).count();
        });
      }
      tg.wait();

      // Pool all cuts from all separators.
      struct PooledCut {
        size_t sep_idx;
        size_t cut_idx;
        double violation;
      };
      std::vector<PooledCut> pool;
      for (size_t i = 0; i < separators_.size(); ++i) {
        const std::string sep_name = separators_[i]->name();
        auto& stats = separator_stats_[sep_name];
        stats.time_seconds += elapsed[i];
        if (!results[i].empty()) stats.rounds_called++;
        // Per-separator cap (optional pre-filter)
        const auto it_cap = per_separator_max_cuts_.find(sep_name);
        const int32_t sep_cap = (it_cap != per_separator_max_cuts_.end())
                                    ? it_cap->second
                                    : 0;  // 0 = unlimited per separator

        auto& cuts = results[i];
        std::stable_sort(cuts.begin(), cuts.end(),
                         [](const sep::Cut& a, const sep::Cut& b) {
                           return a.violation > b.violation;
                         });
        int32_t sep_added = 0;
        for (size_t j = 0; j < cuts.size(); ++j) {
          if (sep_cap > 0 && sep_added >= sep_cap) break;
          pool.push_back({i, j, cuts[j].violation});
          sep_added++;
        }
      }

      // Sort entire pool by violation descending, add top max_cuts_per_round_
      std::stable_sort(pool.begin(), pool.end(),
                       [](const PooledCut& a, const PooledCut& b) {
                         return a.violation > b.violation;
                       });
      int32_t round_added = 0;
      for (const auto& pc : pool) {
        if (max_cuts_per_round_ > 0 && round_added >= max_cuts_per_round_)
          break;
        auto& cut = results[pc.sep_idx][pc.cut_idx];
        const std::string& sep_name = separators_[pc.sep_idx]->name();
        auto& stats = separator_stats_[sep_name];
        std::vector<HighsInt> hi(cut.indices.begin(), cut.indices.end());
        // Disable domain propagation for user cuts: fractional
        // coefficients (e.g. RCI's 2*q_i/Q_r) can cause HiGHS
        // propagation to over-tighten bounds and prune the
        // optimal solution.
        cutpool.addCut(mipsolver, hi.data(), cut.values.data(),
                       static_cast<HighsInt>(hi.size()), cut.rhs,
                       /*integral=*/false,
                       /*propagate=*/false);
        stats.cuts_added++;
        total_cuts_++;
        round_added++;
      }
    }
    separation_rounds_++;
  });

  // Register feasibility check for incumbent solutions.
  // HiGHS heuristics (feasibility pump, rounding, sub-MIPs) may produce
  // integer solutions that bypass the separation loop. This callback
  // rejects any solution with violated SEC/GSEC constraints.
  HighsUserSeparator::setFeasibilityCheck(
      [this](const std::vector<double>& sol) -> bool {
        const int32_t m = num_edges_;
        const int32_t n = num_nodes_;
        if (static_cast<int32_t>(sol.size()) < m + n) return true;

        // Always enforce SEC feasibility on integer incumbents, even when
        // SEC cut generation is disabled, to preserve route connectivity.
        sep::SeparationOracle oracle(prob_);
        return oracle.is_feasible(
            std::span(sol.data(), static_cast<size_t>(m)),
            std::span(sol.data() + m, static_cast<size_t>(n)), x_offset(),
            y_offset());
      });
}

void HiGHSBridge::set_labeling_bounds(std::vector<double> f,
                                      std::vector<double> b,
                                      double correction) {
  fwd_bounds_ = std::move(f);
  bwd_bounds_ = std::move(b);
  correction_ = correction;
}

void HiGHSBridge::set_all_pairs_bounds(std::vector<double> dist) {
  all_pairs_ = std::move(dist);
}

void HiGHSBridge::install_propagator() {
  // No early return: even without bounds, RC fixing (Trigger C/D) is
  // independent.
  const bool have_bounds = !fwd_bounds_.empty() && !bwd_bounds_.empty();
  const bool bounds_prop = bounds_propagation_ && have_bounds;

  const auto& graph = prob_.graph();
  const int32_t m = num_edges_;
  const int32_t n = num_nodes_;
  constexpr double inf = std::numeric_limits<double>::infinity();

  // Shared state for the propagator callback (captured by the lambda).
  // These are mutable across callback invocations.
  auto last_ub = std::make_shared<double>(inf);
  auto processed_edge = std::make_shared<std::vector<bool>>(m, false);
  propagator_fixings_ = std::make_shared<int64_t>(0);
  sweep_fixings_ = std::make_shared<int64_t>(0);
  chain_fixings_ = std::make_shared<int64_t>(0);
  sweep_node_fixings_ = std::make_shared<int64_t>(0);
  chain_node_fixings_ = std::make_shared<int64_t>(0);
  ub_improvements_ = std::make_shared<int64_t>(0);
  propagator_calls_ = std::make_shared<int64_t>(0);
  propagator_time_seconds_ = std::make_shared<double>(0.0);
  rc_fix0_count_ = std::make_shared<int64_t>(0);
  rc_fix1_count_ = std::make_shared<int64_t>(0);
  rc_label_runs_ = std::make_shared<int64_t>(0);
  rc_callback_runs_ = std::make_shared<int64_t>(0);
  rc_time_seconds_ = std::make_shared<double>(0.0);
  auto propagator_fixings = propagator_fixings_;
  auto sweep_fixings = sweep_fixings_;
  auto chain_fixings = chain_fixings_;
  auto sweep_node_fixings = sweep_node_fixings_;
  auto chain_node_fixings = chain_node_fixings_;
  auto ub_improvements = ub_improvements_;
  auto propagator_calls = propagator_calls_;
  auto propagator_time = propagator_time_seconds_;
  auto rc_fix0 = rc_fix0_count_;
  auto rc_fix1 = rc_fix1_count_;
  auto rc_runs = rc_label_runs_;
  auto rc_callbacks = rc_callback_runs_;
  auto rc_time = rc_time_seconds_;
  auto rc_settings = rc_settings_;
  bool obj_is_integer = obj_is_integer_;
  auto rc_adaptive_disabled = std::make_shared<bool>(false);
  auto rc_low_yield_streak = std::make_shared<int32_t>(0);

  // All-pairs bounds for stronger Trigger B (empty if not computed)
  auto ap = std::make_shared<std::vector<double>>(std::move(all_pairs_));

  // Pre-build adjacency: for each node, list of (edge_idx, neighbor).
  auto adj = std::make_shared<std::vector<std::vector<PropagatorAdjEntry>>>(n);
  for (auto e : graph.edges()) {
    int32_t u = graph.edge_source(e);
    int32_t v = graph.edge_target(e);
    (*adj)[u].push_back({e, v});
    (*adj)[v].push_back({e, u});
  }

  // Pre-build edge costs vector for fast access
  auto edge_costs = std::make_shared<std::vector<double>>(m);
  for (auto e : graph.edges()) {
    (*edge_costs)[e] = prob_.edge_cost(e);
  }

  // Pre-build profit vector
  auto profits = std::make_shared<std::vector<double>>(n);
  for (int32_t i = 0; i < n; ++i) {
    (*profits)[i] = prob_.profit(i);
  }

  // Capture bounds by value for the propagator lambda (immutable after init)
  auto prop_fwd = fwd_bounds_;
  auto prop_bwd = bwd_bounds_;
  auto prop_corr = correction_;

  HighsUserPropagator::setCallback([=, &prob = prob_](
                                       HighsDomain& domain,
                                       const HighsMipSolver& mipsolver,
                                       const HighsLpRelaxation& lp) {
    // Skip sub-MIPs: different column space
    if (mipsolver.submip) return;

    const auto& f = prop_fwd;
    const auto& b = prop_bwd;
    const double corr = prop_corr;

    double ub = mipsolver.mipdata_->upper_limit;
    if (ub >= inf) return;

    (*propagator_calls)++;
    const auto propagator_t0 = std::chrono::steady_clock::now();
    const auto& ec = *edge_costs;
    const auto& pr = *profits;
    const auto& adjacency = *adj;
    auto& proc = *processed_edge;
    int64_t& fixings = *propagator_fixings;
    auto fix_nodes_with_no_open_edges =
        [&](bool skip_terminals, bool count_toward_total_fixings,
            const std::shared_ptr<int64_t>& node_fix_counter) {
          for (int32_t i = 0; i < n; ++i) {
            if (domain.col_upper_[m + i] < 0.5) continue;
            if (skip_terminals && (i == prob.source() || i == prob.target())) {
              continue;
            }
            bool all_fixed = true;
            for (const auto& [e, nb] : adjacency[i]) {
              if (domain.col_upper_[e] > 0.5) {
                all_fixed = false;
                break;
              }
            }
            if (!all_fixed) continue;
            domain.changeBound(HighsBoundType::kUpper, m + i, 0.0,
                               HighsDomain::Reason::unspecified());
            if (count_toward_total_fixings) {
              fixings++;
            }
            (*node_fix_counter)++;
          }
        };

    bool ub_improved = (ub < *last_ub - 1e-9);
    if (ub_improved) {
      *last_ub = ub;
      (*ub_improvements)++;
    }

    // ── Trigger A: UB improved → full sweep ──
    // (requires bounds-based propagation)
    if (bounds_prop && ub_improved) {
      // Reset processed edges — tighter UB may enable new fixings
      std::fill(proc.begin(), proc.end(), false);

      for (int32_t e = 0; e < m; ++e) {
        if (domain.col_upper_[e] < 0.5) continue;  // already fixed

        // Try both orientations
        int32_t u = prob.graph().edge_source(e);
        int32_t v = prob.graph().edge_target(e);

        double lb1 = f[u] + ec[e] + b[v] + corr;
        double lb2 = f[v] + ec[e] + b[u] + corr;
        double lb = std::min(lb1, lb2);

        if (lb > ub + 1e-6) {
          domain.changeBound(HighsBoundType::kUpper, e, 0.0,
                             HighsDomain::Reason::unspecified());
          fixings++;
          (*sweep_fixings)++;
        }
      }

      fix_nodes_with_no_open_edges(false, false, sweep_node_fixings);
    }

    // ── Trigger B: Fixed edge x_(a,i) = 1 → chained bounds ──
    // (requires bounds-based propagation)
    if (bounds_prop) {
      const auto& apd = *ap;
      const bool have_all_pairs = !apd.empty();
      const int32_t depot = prob.depot();
      const int64_t chain_before = *chain_fixings;

      for (int32_t e = 0; e < m; ++e) {
        if (proc[e]) continue;
        if (domain.col_lower_[e] < 0.5) continue;  // not fixed to 1
        proc[e] = true;

        int32_t a = prob.graph().edge_source(e);
        int32_t i = prob.graph().edge_target(e);

        if (have_all_pairs) {
          // All-pairs Trigger B: scan ALL unfixed edges
          const double d_depot_a = apd[depot * n + a];
          const double d_depot_i = apd[depot * n + i];
          const double d_i_depot = apd[i * n + depot];
          const double d_a_depot = apd[a * n + depot];
          const double cost_ai = ec[e];

          for (int32_t ej = 0; ej < m; ++ej) {
            if (ej == e) continue;
            if (domain.col_upper_[ej] < 0.5) continue;

            int32_t u = prob.graph().edge_source(ej);
            int32_t v = prob.graph().edge_target(ej);

            // Orientation 1: depot→...→a→i→...→u→v→...→depot
            double lb1 = d_depot_a + cost_ai + apd[i * n + u] + ec[ej] +
                         apd[v * n + depot] + corr;
            // Orientation 2: depot→...→u→v→...→a→i→...→depot
            double lb2 = apd[depot * n + u] + ec[ej] + apd[v * n + a] +
                         cost_ai + d_i_depot + corr;
            // Also try reversed edge orientation (u,v) as (v,u)
            double lb3 = d_depot_a + cost_ai + apd[i * n + v] + ec[ej] +
                         apd[u * n + depot] + corr;
            double lb4 = apd[depot * n + v] + ec[ej] + apd[u * n + a] +
                         cost_ai + d_i_depot + corr;
            // And reversed fixed edge: depot→...→i→a→...
            double lb5 = d_depot_i + cost_ai + apd[a * n + u] + ec[ej] +
                         apd[v * n + depot] + corr;
            double lb6 = apd[depot * n + u] + ec[ej] + apd[v * n + i] +
                         cost_ai + d_a_depot + corr;
            double lb7 = d_depot_i + cost_ai + apd[a * n + v] + ec[ej] +
                         apd[u * n + depot] + corr;
            double lb8 = apd[depot * n + v] + ec[ej] + apd[u * n + i] +
                         cost_ai + d_a_depot + corr;

            double lb = std::min({lb1, lb2, lb3, lb4, lb5, lb6, lb7, lb8});

            if (lb > ub + 1e-6) {
              domain.changeBound(HighsBoundType::kUpper, ej, 0.0,
                                 HighsDomain::Reason::unspecified());
              fixings++;
              (*chain_fixings)++;
            }
          }
        } else {
          // Fallback: neighbor-only scan (original Trigger B)
          double cost_a_to_i = f[a] + ec[e] - pr[i];

          for (const auto& [ej, j] : adjacency[i]) {
            if (ej == e) continue;
            if (domain.col_upper_[ej] < 0.5) continue;

            double lb = cost_a_to_i + ec[ej] + b[j] + corr;
            if (lb > ub + 1e-6) {
              domain.changeBound(HighsBoundType::kUpper, ej, 0.0,
                                 HighsDomain::Reason::unspecified());
              fixings++;
              (*chain_fixings)++;
            }
          }

          double cost_via_i_return = ec[e] + b[i] - pr[a];

          for (const auto& [ek, k] : adjacency[a]) {
            if (ek == e) continue;
            if (domain.col_upper_[ek] < 0.5) continue;

            double lb = f[k] + ec[ek] + cost_via_i_return + corr;
            if (lb > ub + 1e-6) {
              domain.changeBound(HighsBoundType::kUpper, ek, 0.0,
                                 HighsDomain::Reason::unspecified());
              fixings++;
              (*chain_fixings)++;
            }
          }
        }
      }

      // Fix nodes where all incident edges are now fixed to 0.
      if (*chain_fixings > chain_before) {
        fix_nodes_with_no_open_edges(true, true, chain_node_fixings);
      }
    }  // if (bounds_prop)

    // ── Trigger C: Lagrangian reduced-cost fixing (edges → 0) ──
    // ── Trigger D: Lagrangian reduced-cost fixing (nodes → 1) ──
    if (rc_settings.strategy != RCFixingStrategy::off) {
      const bool run_rc =
          should_run_rc_fixing(rc_settings, ub_improved, *ub_improvements,
                               *propagator_calls, *rc_adaptive_disabled);

      if (run_rc && lp.getStatus() == HighsLpRelaxation::Status::kOptimal) {
        auto rc_t0 = std::chrono::steady_clock::now();
        int64_t fix0_before = *rc_fix0;
        int64_t fix1_before = *rc_fix1;
        const auto& lp_sol = lp.getSolution();
        double z_LP = lp.getObjective();
        if (obj_is_integer) {
          z_LP = std::ceil(z_LP - 1e-9);
        }
        int32_t ncols = static_cast<int32_t>(lp_sol.col_dual.size());

        // Build reduced-cost edge weights and node profits.
        // RC for variable j = c_j - dual info. For MIP columns at LP:
        //   rc_j = col_dual[j] (HiGHS stores reduced costs in col_dual)
        // The Lagrangian edge cost = original_cost + rc (reduced cost acts as
        // perturbation). Actually, the reduced cost already IS the Lagrangian
        // cost relative to the LP basis. For labeling: use rc_j directly as the
        // edge weight. For y-vars: rc_{m+i} as the negative profit
        // contribution.
        if (ncols >= m + n) {
          std::vector<double> rc_edge_costs(m);
          std::vector<double> rc_profits(n);
          for (int32_t e = 0; e < m; ++e) {
            rc_edge_costs[e] = lp_sol.col_dual[e];
          }
          for (int32_t i = 0; i < n; ++i) {
            // y_i has objective coefficient -profit(i), so
            // rc for y_i = col_dual[m+i]. The labeling subtracts profits,
            // so rc_profit = -col_dual[m+i] to match the labeling convention.
            rc_profits[i] = -lp_sol.col_dual[m + i];
          }

          // Correction under RC costs: for tours, depot RC profit
          double rc_corr = prob.is_tour() ? rc_profits[prob.source()] : 0.0;

          // Forward + backward labeling with RC costs (parallel)
          std::vector<double> rc_fwd, rc_bwd;
          {
            tbb::task_group tg;
            tg.run([&] {
              rc_fwd = preprocess::labeling_from(prob, prob.source(),
                                                 rc_edge_costs, rc_profits);
            });
            if (prob.is_tour()) {
              tg.wait();
              rc_bwd = rc_fwd;
            } else {
              tg.run([&] {
                rc_bwd = preprocess::labeling_from(prob, prob.target(),
                                                   rc_edge_costs, rc_profits);
              });
              tg.wait();
            }
          }
          (*rc_runs) += prob.is_tour() ? 1 : 2;

          // Compute z_LR: cheapest resource-feasible tour/path under RC costs
          double z_LR = inf;
          if (prob.is_tour()) {
            // Tour: z_LR = min over edges (rc_fwd[u] + rc_e + rc_fwd[v]) +
            // rc_corr
            for (int32_t e = 0; e < m; ++e) {
              int32_t u = prob.graph().edge_source(e);
              int32_t v = prob.graph().edge_target(e);
              if (rc_fwd[u] < inf && rc_fwd[v] < inf) {
                double val = rc_fwd[u] + rc_edge_costs[e] + rc_fwd[v] + rc_corr;
                z_LR = std::min(z_LR, val);
              }
            }
          } else {
            // Path: z_LR = rc_fwd[target]
            z_LR = rc_fwd[prob.target()];
          }

          // Guard: skip RC fixing when the Lagrangian gap
          // (z_LP − z_LR) is negative or exceeds the
          // optimality gap.  A negative gap means z_LR > z_LP
          // (labeling looser than LP); a large positive gap
          // means the non-elementary excess overestimates
          // the true per-edge penalty.
          double lagrangian_gap = z_LP - z_LR;
          double optimality_gap = ub - z_LP;
          if (z_LR < inf && lagrangian_gap >= 0.0 &&
              lagrangian_gap < optimality_gap + 1e-6) {
            // Trigger C: Fix edges to 0
            for (int32_t e = 0; e < m; ++e) {
              if (domain.col_upper_[e] < 0.5 || domain.col_lower_[e] > 0.5)
                continue;

              int32_t u = prob.graph().edge_source(e);
              int32_t v = prob.graph().edge_target(e);

              double lb1 =
                  (rc_fwd[u] < inf && rc_bwd[v] < inf)
                      ? rc_fwd[u] + rc_edge_costs[e] + rc_bwd[v] + rc_corr
                      : inf;
              double lb2 =
                  (rc_fwd[v] < inf && rc_bwd[u] < inf)
                      ? rc_fwd[v] + rc_edge_costs[e] + rc_bwd[u] + rc_corr
                      : inf;
              double edge_lb = std::min(lb1, lb2);

              double excess = edge_lb - z_LR;
              if (z_LP + excess > ub + 1e-6) {
                domain.changeBound(HighsBoundType::kUpper, e, 0.0,
                                   HighsDomain::Reason::unspecified());
                (*rc_fix0)++;
              }
            }

            // Fix nodes where all incident edges are now fixed to 0
            for (int32_t i = 0; i < n; ++i) {
              if (domain.col_upper_[m + i] < 0.5) continue;
              if (i == prob.source() || i == prob.target()) continue;
              bool all_fixed = true;
              for (const auto& [e, nb] : adjacency[i]) {
                if (domain.col_upper_[e] > 0.5) {
                  all_fixed = false;
                  break;
                }
              }
              if (all_fixed) {
                domain.changeBound(HighsBoundType::kUpper, m + i, 0.0,
                                   HighsDomain::Reason::unspecified());
              }
            }

            // Trigger D: Fix nodes to 1 (optional, expensive)
            if (rc_settings.fix_to_one) {
              // For each unfixed node, run labeling with that node forbidden
              std::vector<int32_t> candidates;
              for (int32_t i = 0; i < n; ++i) {
                if (i == prob.source() || i == prob.target()) continue;
                if (domain.col_lower_[m + i] > 0.5)
                  continue;  // already fixed to 1
                if (domain.col_upper_[m + i] < 0.5)
                  continue;  // already fixed to 0
                candidates.push_back(i);
              }

              // Parallel labeling with each node forbidden
              std::vector<double> z_LR_minus(candidates.size(), inf);
              tbb::parallel_for(
                  static_cast<size_t>(0), candidates.size(), [&](size_t idx) {
                    int32_t node_i = candidates[idx];
                    // Copy profits and set forbidden node to -1e30
                    std::vector<double> mod_profits = rc_profits;
                    mod_profits[node_i] = -1e30;

                    auto fwd_i = preprocess::labeling_from(
                        prob, prob.source(), rc_edge_costs, mod_profits);

                    if (prob.is_tour()) {
                      double best = inf;
                      for (int32_t e = 0; e < m; ++e) {
                        int32_t u = prob.graph().edge_source(e);
                        int32_t v = prob.graph().edge_target(e);
                        if (fwd_i[u] < inf && fwd_i[v] < inf) {
                          double val =
                              fwd_i[u] + rc_edge_costs[e] + fwd_i[v] + rc_corr;
                          best = std::min(best, val);
                        }
                      }
                      z_LR_minus[idx] = best;
                    } else {
                      z_LR_minus[idx] = fwd_i[prob.target()];
                    }
                  });

              (*rc_runs) += static_cast<int64_t>(candidates.size());

              for (size_t idx = 0; idx < candidates.size(); ++idx) {
                if (z_LR_minus[idx] >= inf) continue;
                double gap = z_LR_minus[idx] - z_LR;
                if (z_LP + gap > ub + 1e-6) {
                  int32_t node_i = candidates[idx];
                  domain.changeBound(HighsBoundType::kLower, m + node_i, 1.0,
                                     HighsDomain::Reason::unspecified());
                  (*rc_fix1)++;
                }
              }
            }
          }
        }

        (*rc_callbacks)++;
        auto rc_t1 = std::chrono::steady_clock::now();
        *rc_time += std::chrono::duration<double>(rc_t1 - rc_t0).count();

        if (rc_settings.strategy == RCFixingStrategy::adaptive &&
            *rc_callbacks >= 2) {
          const int64_t delta_fix =
              (*rc_fix0 - fix0_before) + (*rc_fix1 - fix1_before);
          if (delta_fix < 5) {
            (*rc_low_yield_streak)++;
          } else {
            *rc_low_yield_streak = 0;
          }

          // Two low-yield UB-improvement rounds in a row => disable.
          if (*rc_low_yield_streak >= 2) {
            *rc_adaptive_disabled = true;
          }
        }
      }
    }
    if (propagator_time) {
      const auto propagator_t1 = std::chrono::steady_clock::now();
      *propagator_time +=
          std::chrono::duration<double>(propagator_t1 - propagator_t0).count();
    }
  });
}

void HiGHSBridge::install_heuristic_callback() {
  const bool heuristic_enabled = heuristic_callback_;
  const bool interrupt_enabled = static_cast<bool>(interrupt_flag_);
  const bool logging_enabled = logger_.enabled();
  if (!heuristic_enabled && !interrupt_enabled && !logging_enabled) return;

  auto calls = heuristic_calls_;
  auto solutions = heuristic_solutions_;
  auto heuristic_time_seconds = heuristic_time_seconds_;
  int strategy = heuristic_strategy_;
  int64_t node_interval = std::max<int64_t>(1, heuristic_node_interval_);
  int deterministic_restarts =
      std::max<int32_t>(1, heuristic_deterministic_restarts_);
  // Rate-limit: run heuristic at most once per N nodes.
  auto last_node_count = std::make_shared<int64_t>(0);

  auto cached_x = cached_x_lp_;
  auto cached_y = cached_y_lp_;
  auto cache_mtx = lp_cache_mutex_;
  double lpg_edge_threshold = heu_lpg_edge_threshold_;
  double lpg_node_threshold = heu_lpg_node_threshold_;
  double lpg_lp_threshold = heu_lpg_lp_threshold_;
  double lpg_seed_threshold = heu_lpg_seed_threshold_;
  highs_.setCallback(
      [this, heuristic_enabled, cached_x, cached_y, cache_mtx, calls, solutions,
       heuristic_time_seconds, last_node_count, strategy, node_interval,
       deterministic_restarts, lpg_edge_threshold, lpg_node_threshold,
       lpg_lp_threshold,
       lpg_seed_threshold](int callback_type, const std::string& message,
                           const HighsCallbackOutput* data_out,
                           HighsCallbackInput* data_in, void* /*user_data*/) {
        if (callback_type ==
            static_cast<int>(HighsCallbackType::kCallbackLogging)) {
          if (!logger_.enabled()) return;
          logger_.log(std::string_view{message});
          return;
        }

        if (callback_type ==
            static_cast<int>(HighsCallbackType::kCallbackMipInterrupt)) {
          if (interrupt_flag_ &&
              interrupt_flag_->load(std::memory_order_relaxed)) {
            data_in->user_interrupt = true;
            return;
          }
          if (!heuristic_enabled) return;

          // Only stop when HiGHS has already accepted an incumbent.
          if (std::isfinite(data_out->mip_primal_bound) &&
              std::isfinite(data_out->mip_dual_bound) &&
              data_out->mip_primal_bound <= data_out->mip_dual_bound + 1e-6) {
            data_in->user_interrupt = true;
          }
          return;
        }

        if (!heuristic_enabled) return;

        // Only handle kCallbackMipUserSolution below.
        if (callback_type !=
            static_cast<int>(HighsCallbackType::kCallbackMipUserSolution)) {
          return;
        }

        // Skip the initial query after setup (pre-solve heuristic already ran)
        if (data_out->external_solution_query_origin ==
            kExternalMipSolutionQueryOriginAfterSetup) {
          return;
        }

        // Rate-limit: skip if not enough nodes since last call
        int64_t current_nodes = data_out->mip_node_count;
        if (current_nodes - *last_node_count < node_interval) {
          return;
        }
        *last_node_count = current_nodes;

        // Read cached LP values from separator callback
        std::vector<double> x_lp, y_lp;
        {
          std::lock_guard lock(*cache_mtx);
          if (cached_x->empty()) return;
          x_lp = *cached_x;
          y_lp = *cached_y;
        }

        // Build incumbent vector from current best solution (if available)
        std::vector<double> incumbent;
        if (!data_out->mip_solution.empty() &&
            data_out->mip_primal_bound < 1e20) {
          incumbent = data_out->mip_solution;
        }

        (*calls)++;
        const auto heuristic_t0 = std::chrono::steady_clock::now();
        heuristic::HeuristicResult result = heuristic::lp_guided_heuristic(
            prob_, x_lp, y_lp, incumbent, 0.0, strategy, deterministic_restarts,
            static_cast<uint32_t>(current_nodes), lpg_edge_threshold,
            lpg_node_threshold, lpg_lp_threshold, lpg_seed_threshold);

        if (!result.col_values.empty() &&
            result.objective < data_out->mip_primal_bound - 1e-6) {
          data_in->user_has_solution = true;
          data_in->user_solution = std::move(result.col_values);
          (*solutions)++;
        }
        if (heuristic_time_seconds) {
          const auto heuristic_t1 = std::chrono::steady_clock::now();
          *heuristic_time_seconds +=
              std::chrono::duration<double>(heuristic_t1 - heuristic_t0)
                  .count();
        }
      },
      nullptr);

  if (heuristic_enabled) {
    highs_.startCallback(HighsCallbackType::kCallbackMipUserSolution);
  }
  highs_.startCallback(HighsCallbackType::kCallbackMipInterrupt);
  if (logging_enabled) {
    highs_.startCallback(HighsCallbackType::kCallbackLogging);
  }
}

std::vector<int32_t> HiGHSBridge::order_path(
    const std::vector<int32_t>& visited_nodes,
    const std::vector<int32_t>& active_edges) const {
  if (visited_nodes.empty()) return {};

  const auto& graph = prob_.graph();
  const int32_t start = prob_.source();

  // Build adjacency from active edges (undirected)
  std::vector<std::vector<int32_t>> adj(num_nodes_);
  for (int32_t edge_idx : active_edges) {
    int32_t u = graph.edge_source(edge_idx);
    int32_t v = graph.edge_target(edge_idx);
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  // Follow the chain from source
  std::vector<int32_t> ordered;
  ordered.reserve(visited_nodes.size() + 1);

  std::vector<bool> seen(num_nodes_, false);
  int32_t current = start;

  do {
    ordered.push_back(current);
    seen[current] = true;
    int32_t next = -1;
    for (int32_t nb : adj[current]) {
      if (!seen[nb]) {
        next = nb;
        break;
      }
    }
    current = next;
  } while (current != -1);

  // For tours, close the loop
  if (prob_.is_tour() && ordered.size() > 1) {
    ordered.push_back(start);
  }

  // If some visited nodes weren't reached (disconnected), append them
  for (int32_t v : visited_nodes) {
    if (!seen[v]) {
      ordered.push_back(v);
    }
  }

  return ordered;
}

SolveResult HiGHSBridge::extract_result() const {
  SolveResult result;

  auto model_status = highs_.getModelStatus();
  switch (model_status) {
    case HighsModelStatus::kNotset:
      // In some callback-heavy runs HiGHS can return kNotset at exit;
      // treat this as an interrupted run and recover incumbent/bounds.
      result.status = SolveResult::Status::TimeLimit;
      break;
    case HighsModelStatus::kOptimal:
      result.status = SolveResult::Status::Optimal;
      break;
    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnboundedOrInfeasible:
      result.status = SolveResult::Status::Infeasible;
      return result;
    case HighsModelStatus::kUnbounded:
      result.status = SolveResult::Status::Unbounded;
      return result;
    case HighsModelStatus::kModelEmpty:
      result.status = SolveResult::Status::Optimal;
      break;
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
      result.status = SolveResult::Status::Feasible;
      break;
    case HighsModelStatus::kInterrupt:
    case HighsModelStatus::kHighsInterrupt:
      result.status = SolveResult::Status::Feasible;
      break;
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kSolutionLimit:
    case HighsModelStatus::kMemoryLimit:
    case HighsModelStatus::kUnknown:
      result.status = SolveResult::Status::TimeLimit;
      break;
    default:
      result.status = SolveResult::Status::Error;
      return result;
  }

  highs_.getInfoValue("objective_function_value", result.objective);

  int64_t mip_nodes = 0;
  if (highs_.getInfoValue("mip_node_count", mip_nodes) == HighsStatus::kOk) {
    result.nodes = mip_nodes;
  }
  HighsInt simplex_iters = -1;
  if (highs_.getInfoValue("simplex_iteration_count", simplex_iters) ==
      HighsStatus::kOk) {
    result.simplex_iterations = static_cast<int64_t>(simplex_iters);
  }
  highs_.getInfoValue("mip_dual_bound", result.bound);
  double mip_gap = 0.0;
  highs_.getInfoValue("mip_gap", mip_gap);
  result.gap = mip_gap / 100.0;

  if (result.status == SolveResult::Status::Feasible &&
      std::isfinite(result.objective) && std::isfinite(result.bound) &&
      result.objective <= result.bound + 1e-6) {
    result.status = SolveResult::Status::Optimal;
    result.gap = 0.0;
  }

  // Extract visited nodes and active edges
  const auto& sol = highs_.getSolution();
  const auto& graph = prob_.graph();

  std::vector<int32_t> visited_nodes;
  for (int32_t i = 0; i < num_nodes_; ++i) {
    if (sol.col_value[num_edges_ + i] > 0.5) {
      visited_nodes.push_back(i);
    }
  }

  std::vector<int32_t> active_edges;
  for (auto e : graph.edges()) {
    if (sol.col_value[e] > 0.5) {
      active_edges.push_back(e);
    }
  }

  // Order the tour by following edges from depot
  result.tour = order_path(visited_nodes, active_edges);
  result.tour_arcs = std::move(active_edges);

  // Print propagator statistics
  if (propagator_calls_) {
    const int64_t sweep_nodes = sweep_node_fixings_ ? *sweep_node_fixings_ : 0;
    const int64_t sweep_total = *sweep_fixings_ + sweep_nodes;
    const int64_t chain_nodes = chain_node_fixings_ ? *chain_node_fixings_ : 0;
    const int64_t chain_total =
        (chain_fixings_ ? *chain_fixings_ : 0) + chain_nodes;
    const int64_t propagator_fixings = sweep_total + chain_total;
    logger_.log(
        "Propagator: calls={}  fixings={} (sweep={} [edges={} nodes={}] "
        "chain={} [edges={} nodes={}])  time={:.3f}s",
        *propagator_calls_, propagator_fixings, sweep_total, *sweep_fixings_,
        sweep_nodes, chain_total, *chain_fixings_, chain_nodes,
        propagator_time_seconds_ ? *propagator_time_seconds_ : 0.0);
  }
  if (rc_fix0_count_ && (*rc_fix0_count_ > 0 || *rc_fix1_count_ > 0)) {
    const int64_t rc_fixings = *rc_fix0_count_ + *rc_fix1_count_;
    const int64_t rc_calls = rc_callback_runs_ ? *rc_callback_runs_ : 0;
    logger_.log(
        "RC fixing:  calls={}  fixings={} (edges={} nodes={})  "
        "labeling_runs={}  time={:.3f}s",
        rc_calls, rc_fixings, *rc_fix0_count_, *rc_fix1_count_, *rc_label_runs_,
        rc_time_seconds_ ? *rc_time_seconds_ : 0.0);
  }

  // Print heuristic callback statistics for evaluation.
  const int64_t heuristic_calls = heuristic_calls_ ? *heuristic_calls_ : 0;
  const int64_t heuristic_solutions =
      heuristic_solutions_ ? *heuristic_solutions_ : 0;
  const double heuristic_time =
      heuristic_time_seconds_ ? *heuristic_time_seconds_ : 0.0;
  if (heuristic_callback_ || heuristic_calls > 0) {
    logger_.log("Heuristic:  calls={}  solutions={}  time={:.3f}s",
                heuristic_calls, heuristic_solutions, heuristic_time);
  }

  // Attach separator statistics
  result.separator_stats = separator_stats_;
  result.total_cuts = total_cuts_;
  result.separation_rounds = separation_rounds_;

  return result;
}

}  // namespace cptp
