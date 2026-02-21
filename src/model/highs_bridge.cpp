#include "model/highs_bridge.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <mutex>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "mip/HighsUserSeparator.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/reachability.h"
#include "sep/sec_separator.h"
#include "sep/separation_context.h"

#include <tbb/task_group.h>

namespace cptp {

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
    if (upper_bound_ < std::numeric_limits<double>::infinity()) {
        auto eliminated = preprocess::edge_elimination(prob_, upper_bound_);
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
            // Skip sub-MIPs: they have a different column space
            if (mipsolver.submip) return;

            // Amortized separation: skip rounds to reduce overhead
            int64_t call = separation_calls_++;
            if (separation_interval_ > 1 && (call % separation_interval_) != 0)
                return;

            const auto& sol = lpRelaxation.getSolution();
            const int32_t m = num_edges_;
            const int32_t n = num_nodes_;

            std::vector<double> x_vals(sol.col_value.begin(),
                                        sol.col_value.begin() + m);
            std::vector<double> y_vals(sol.col_value.begin() + m,
                                        sol.col_value.begin() + m + n);

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

            // Run all separators in parallel, timing each one
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

                // Sort by violation (descending) and limit count
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

    // Attach separator statistics
    result.separator_stats = separator_stats_;
    result.total_cuts = total_cuts_;
    result.separation_rounds = separation_rounds_;

    return result;
}

}  // namespace cptp
