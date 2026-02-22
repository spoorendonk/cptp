#include "model/highs_bridge.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <mutex>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsUserPropagator.h"
#include "mip/HighsUserSeparator.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/reachability.h"
#include "sep/sec_separator.h"
#include "sep/separation_context.h"

#include <tbb/task_group.h>

namespace rcspp {

HiGHSBridge::HiGHSBridge(const Problem& prob, Highs& highs, double frac_tol)
    : prob_(prob),
      highs_(highs),
      num_edges_(prob.num_edges()),
      num_nodes_(prob.num_nodes()),
      frac_tol_(frac_tol),
      int_tol_(1e-6) {
    // Integral feasibility: tight but not zero to avoid rejecting
    // valid solutions due to floating-point noise in LP values.
    // Fractional separation uses a looser tolerance (default 1e-2).
}

HiGHSBridge::~HiGHSBridge() {
    HighsUserSeparator::clearCallback();
    HighsUserPropagator::clearCallback();
}

void HiGHSBridge::add_separator(std::unique_ptr<sep::Separator> sep) {
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

    // Fix source and target y variables to 1 via bounds
    col_lower[m + prob_.source()] = 1.0;
    col_lower[m + prob_.target()] = 1.0;

    // Fix unreachable nodes: demand-reachability preprocessing
    if (prob_.capacity() > 0 && prob_.capacity() < 1e17) {
        auto reachable = preprocess::demand_reachability(prob_);
        for (int32_t i = 0; i < n; ++i) {
            if (!reachable[i] && i != prob_.source() && i != prob_.target()) {
                col_upper[m + i] = 0.0;
                for (auto e : graph.incident_edges(i))
                    col_upper[e] = 0.0;
            }
        }
    }

    // Edge elimination: capacity-aware 2-cycle labeling bounds
    if (upper_bound_ < std::numeric_limits<double>::infinity()
        && !fwd_bounds_.empty() && !bwd_bounds_.empty()) {
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
        for (int32_t i = 0; i < n; ++i) {
            if (i == prob_.source() || i == prob_.target() || col_upper[m + i] == 0.0) continue;
            bool all_eliminated = true;
            for (auto e : graph.incident_edges(i)) {
                if (col_upper[e] > 0.0) { all_eliminated = false; break; }
            }
            if (all_eliminated) {
                col_upper[m + i] = 0.0;
                node_elim_count++;
            }
        }
        if (edge_elim_count > 0 || node_elim_count > 0) {
            std::cerr << "Edge elimination: " << edge_elim_count << "/" << m
                      << " edges, " << node_elim_count << "/" << (n - 1)
                      << " nodes fixed to 0\n";
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

    // Add columns with costs
    for (int i = 0; i < total_vars; ++i) {
        highs_.addCol(col_cost[i], col_lower[i], col_upper[i],
                       0, nullptr, nullptr);
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
        bool is_terminal = !prob_.is_tour() &&
                           (i == prob_.source() || i == prob_.target());
        indices.push_back(m + i);
        values.push_back(is_terminal ? -1.0 : -2.0);
        highs_.addRow(0.0, 0.0, static_cast<int>(indices.size()),
                       indices.data(), values.data());
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
                           static_cast<int>(indices.size()),
                           indices.data(), values.data());
        }
    }
}

void HiGHSBridge::install_separators() {
    // Capture what the callback needs by pointer/reference.
    // The bridge must outlive the highs.run() call.
    HighsUserSeparator::setCallback(
        [this](const HighsLpRelaxation& lpRelaxation,
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
            };

            const bool is_submip = !orig_to_reduced.empty();

            if (is_submip) {
                // Sub-MIP: only run SEC (lazy constraint) to keep overhead low.
                // Other separators (RCI, Multistar, etc.) add tightening cuts
                // but aren't needed for feasibility in short-lived sub-MIPs.
                sep::SECSeparator sec;
                auto cuts = sec.separate(ctx);

                std::sort(cuts.begin(), cuts.end(),
                    [](const sep::Cut& a, const sep::Cut& b) {
                        return a.violation > b.violation;
                    });
                int32_t limit = (max_cuts_per_sep_ > 0)
                    ? std::min(static_cast<int32_t>(cuts.size()), max_cuts_per_sep_)
                    : static_cast<int32_t>(cuts.size());

                for (int32_t j = 0; j < limit; ++j) {
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
                            double fixed_val = (orig_col < m)
                                ? x_vals[orig_col]
                                : y_vals[orig_col - m];
                            rhs -= cut.values[k] * fixed_val;
                        }
                    }

                    if (!red_idx.empty()) {
                        cutpool.addCut(mipsolver,
                                       red_idx.data(), red_val.data(),
                                       static_cast<HighsInt>(red_idx.size()),
                                       rhs);
                        total_cuts_++;
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
                        elapsed[i] = std::chrono::duration<double>(
                            clock::now() - t0).count();
                    });
                }
                tg.wait();

                // Add cuts to the cutpool (limited per separator)
                for (size_t i = 0; i < separators_.size(); ++i) {
                    auto& stats = separator_stats_[separators_[i]->name()];
                    stats.time_seconds += elapsed[i];
                    if (!results[i].empty()) stats.rounds_called++;

                    auto& cuts = results[i];
                    std::sort(cuts.begin(), cuts.end(),
                        [](const sep::Cut& a, const sep::Cut& b) {
                            return a.violation > b.violation;
                        });
                    int32_t limit = (max_cuts_per_sep_ > 0)
                        ? std::min(static_cast<int32_t>(cuts.size()), max_cuts_per_sep_)
                        : static_cast<int32_t>(cuts.size());

                    for (int32_t j = 0; j < limit; ++j) {
                        auto& cut = cuts[j];
                        std::vector<HighsInt> hi(cut.indices.begin(), cut.indices.end());
                        cutpool.addCut(mipsolver,
                                       hi.data(), cut.values.data(),
                                       static_cast<HighsInt>(hi.size()),
                                       cut.rhs);
                        stats.cuts_added++;
                        total_cuts_++;
                    }
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

            const auto& graph = prob_.graph();
            const double graph_tol = 1e-6;  // for edge inclusion in support graph

            // Build support graph and Gomory-Hu tree for feasibility check
            digraph_builder builder(n);
            for (auto e : graph.edges()) {
                double xval = sol[e];
                if (xval > graph_tol) {
                    int32_t u = graph.edge_source(e);
                    int32_t v = graph.edge_target(e);
                    builder.add_arc(u, v, xval);
                    builder.add_arc(v, u, xval);
                }
            }
            auto [support_graph, capacity] = builder.build();
            gomory_hu_tree flow_tree(support_graph, capacity, prob_.source());

            sep::SeparationContext ctx{
                .problem = prob_,
                .x_values = std::span(sol.data(), static_cast<size_t>(m)),
                .y_values = std::span(sol.data() + m, static_cast<size_t>(n)),
                .x_offset = x_offset(),
                .y_offset = y_offset(),
                .tol = int_tol_,
                .flow_tree = &flow_tree,
            };

            // Only check SEC — these are the lazy constraints
            sep::SECSeparator sec;
            auto cuts = sec.separate(ctx);
            return cuts.empty();
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
    if (fwd_bounds_.empty() || bwd_bounds_.empty()) return;

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
    ub_improvements_ = std::make_shared<int64_t>(0);
    propagator_calls_ = std::make_shared<int64_t>(0);
    auto propagator_fixings = propagator_fixings_;
    auto sweep_fixings = sweep_fixings_;
    auto chain_fixings = chain_fixings_;
    auto ub_improvements = ub_improvements_;
    auto propagator_calls = propagator_calls_;

    // Capture labeling bounds by value (they won't change during solve).
    auto fwd = std::make_shared<std::vector<double>>(fwd_bounds_);
    auto bwd = std::make_shared<std::vector<double>>(bwd_bounds_);
    double corr = correction_;

    // All-pairs bounds for stronger Trigger B (empty if not computed)
    auto ap = std::make_shared<std::vector<double>>(std::move(all_pairs_));

    // Pre-build adjacency: for each node, list of (edge_idx, neighbor)
    struct AdjEntry { int32_t edge; int32_t neighbor; };
    auto adj = std::make_shared<std::vector<std::vector<AdjEntry>>>(n);
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

    HighsUserPropagator::setCallback(
        [=, &prob = prob_](HighsDomain& domain,
                           const HighsMipSolver& mipsolver,
                           const HighsLpRelaxation& /*lp*/) {
            // Skip sub-MIPs: different column space
            if (mipsolver.submip) return;

            double ub = mipsolver.mipdata_->upper_limit;
            if (ub >= inf) return;

            (*propagator_calls)++;

            const auto& f = *fwd;
            const auto& b = *bwd;
            const auto& ec = *edge_costs;
            const auto& pr = *profits;
            const auto& adjacency = *adj;
            auto& proc = *processed_edge;
            int64_t& fixings = *propagator_fixings;

            bool ub_improved = (ub < *last_ub - 1e-9);

            // ── Trigger A: UB improved → full sweep ──
            if (ub_improved) {
                *last_ub = ub;
                (*ub_improvements)++;
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
                        domain.changeBound(
                            HighsBoundType::kUpper, e, 0.0,
                            HighsDomain::Reason::unspecified());
                        fixings++;
                        (*sweep_fixings)++;
                    }
                }

                // Fix nodes where all incident edges are fixed to 0
                for (int32_t i = 0; i < n; ++i) {
                    if (domain.col_upper_[m + i] < 0.5) continue;
                    bool all_fixed = true;
                    for (const auto& [e, nb] : adjacency[i]) {
                        if (domain.col_upper_[e] > 0.5) {
                            all_fixed = false;
                            break;
                        }
                    }
                    if (all_fixed) {
                        domain.changeBound(
                            HighsBoundType::kUpper, m + i, 0.0,
                            HighsDomain::Reason::unspecified());
                    }
                }
            }

            // ── Trigger B: Fixed edge x_(a,i) = 1 → chained bounds ──
            const auto& apd = *ap;
            const bool have_all_pairs = !apd.empty();
            const int32_t depot = prob.depot();

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
                        double lb1 = d_depot_a + cost_ai + apd[i * n + u] + ec[ej] + apd[v * n + depot] + corr;
                        // Orientation 2: depot→...→u→v→...→a→i→...→depot
                        double lb2 = apd[depot * n + u] + ec[ej] + apd[v * n + a] + cost_ai + d_i_depot + corr;
                        // Also try reversed edge orientation (u,v) as (v,u)
                        double lb3 = d_depot_a + cost_ai + apd[i * n + v] + ec[ej] + apd[u * n + depot] + corr;
                        double lb4 = apd[depot * n + v] + ec[ej] + apd[u * n + a] + cost_ai + d_i_depot + corr;
                        // And reversed fixed edge: depot→...→i→a→...
                        double lb5 = d_depot_i + cost_ai + apd[a * n + u] + ec[ej] + apd[v * n + depot] + corr;
                        double lb6 = apd[depot * n + u] + ec[ej] + apd[v * n + i] + cost_ai + d_a_depot + corr;
                        double lb7 = d_depot_i + cost_ai + apd[a * n + v] + ec[ej] + apd[u * n + depot] + corr;
                        double lb8 = apd[depot * n + v] + ec[ej] + apd[u * n + i] + cost_ai + d_a_depot + corr;

                        double lb = std::min({lb1, lb2, lb3, lb4, lb5, lb6, lb7, lb8});

                        if (lb > ub + 1e-6) {
                            domain.changeBound(
                                HighsBoundType::kUpper, ej, 0.0,
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
                            domain.changeBound(
                                HighsBoundType::kUpper, ej, 0.0,
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
                            domain.changeBound(
                                HighsBoundType::kUpper, ek, 0.0,
                                HighsDomain::Reason::unspecified());
                            fixings++;
                            (*chain_fixings)++;
                        }
                    }
                }
            }
        });

    std::cerr << "Installed domain propagator with labeling bounds\n";
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
        case HighsModelStatus::kOptimal:
            result.status = SolveResult::Status::Optimal;
            break;
        case HighsModelStatus::kInfeasible:
            result.status = SolveResult::Status::Infeasible;
            return result;
        case HighsModelStatus::kUnbounded:
            result.status = SolveResult::Status::Unbounded;
            return result;
        case HighsModelStatus::kObjectiveBound:
        case HighsModelStatus::kObjectiveTarget:
            result.status = SolveResult::Status::Feasible;
            break;
        case HighsModelStatus::kTimeLimit:
        case HighsModelStatus::kIterationLimit:
            result.status = SolveResult::Status::TimeLimit;
            break;
        default:
            result.status = SolveResult::Status::Error;
            return result;
    }

    highs_.getInfoValue("objective_function_value", result.objective);

    highs_.getInfoValue("mip_node_count", result.nodes);
    highs_.getInfoValue("mip_dual_bound", result.bound);
    double mip_gap = 0.0;
    highs_.getInfoValue("mip_gap", mip_gap);
    result.gap = mip_gap / 100.0;

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
    if (propagator_calls_ && *propagator_calls_ > 0) {
        std::cerr << "Propagator: " << *propagator_calls_ << " calls, "
                  << *ub_improvements_ << " UB improvements, "
                  << *propagator_fixings_ << " fixings ("
                  << *sweep_fixings_ << " sweep + "
                  << *chain_fixings_ << " chain)\n";
    }

    // Attach separator statistics
    result.separator_stats = separator_stats_;
    result.total_cuts = total_cuts_;
    result.separation_rounds = separation_rounds_;

    return result;
}

}  // namespace rcspp
