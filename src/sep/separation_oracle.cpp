#include "sep/separation_oracle.h"

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "parallel/parallel.h"
#include "sep/comb_separator.h"
#include "sep/largest_side_cache.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/sec_separator.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace cptp::sep {

SeparationOracle::SeparationOracle(const Problem& prob) : prob_(prob) {}

void SeparationOracle::add_separator(std::unique_ptr<Separator> sep) {
    separators_.push_back(std::move(sep));
}

void SeparationOracle::add_default_separators() {
    separators_.push_back(std::make_unique<SECSeparator>());
    separators_.push_back(std::make_unique<RCISeparator>());
    separators_.push_back(std::make_unique<MultistarSeparator>());
    separators_.push_back(std::make_unique<CombSeparator>());
}

std::vector<Cut> SeparationOracle::separate(std::span<const double> x_values,
                                            std::span<const double> y_values, int32_t x_offset,
                                            int32_t y_offset, double tol) const {
    const auto& graph = prob_.graph();
    const int32_t n = prob_.num_nodes();

    // Build support graph from fractional LP values
    constexpr double kGraphTol = 1e-6;
    digraph_builder builder(n);
    for (auto e : graph.edges()) {
        double xval = x_values[e];
        if (xval > kGraphTol) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            builder.add_arc(u, v, xval);
            builder.add_arc(v, u, xval);
        }
    }
    auto [support_graph, capacity] = builder.build();

    // Compute Gomory-Hu tree (n-1 max-flows via Gusfield)
    gomory_hu_tree flow_tree(support_graph, capacity, prob_.source());
    LargestSideCache largest_sides(support_graph, capacity, prob_.source());

    SeparationContext ctx{
        .problem = prob_,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = x_offset,
        .y_offset = y_offset,
        .tol = tol,
        .flow_tree = &flow_tree,
        .support_graph = &support_graph,
        .support_capacity = &capacity,
        .largest_sides = &largest_sides,
    };

    // Run all separators in parallel
    std::vector<std::vector<Cut>> results(separators_.size());
    parallel::task_group tg;
    for (size_t i = 0; i < separators_.size(); ++i) {
        tg.run([&, i] { results[i] = separators_[i]->separate(ctx); });
    }
    tg.wait();

    // Collect, sort by violation, limit per separator
    std::vector<Cut> all_cuts;
    for (size_t i = 0; i < separators_.size(); ++i) {
        auto& cuts = results[i];
        std::ranges::sort(cuts,
                          [](const Cut& a, const Cut& b) { return a.violation > b.violation; });
        int32_t limit = (max_cuts_per_sep_ > 0)
                            ? std::min(static_cast<int32_t>(cuts.size()), max_cuts_per_sep_)
                            : static_cast<int32_t>(cuts.size());
        for (int32_t j = 0; j < limit; ++j) {
            all_cuts.push_back(std::move(cuts[j]));
        }
    }

    // Sort combined results by violation
    std::ranges::sort(all_cuts,
                      [](const Cut& a, const Cut& b) { return a.violation > b.violation; });

    return all_cuts;
}

namespace {

/// The incumbent check's tolerance: is_feasible_by_max_flow's support-graph
/// threshold and SEC violation tolerance, and integral_connectivity()'s
/// integrality tolerance.
constexpr double kIncumbentTol = 1e-6;

}  // namespace

bool SeparationOracle::is_feasible(std::span<const double> x_values,
                                   std::span<const double> y_values, int32_t x_offset,
                                   int32_t y_offset) const {
    if (auto connected = integral_connectivity(x_values, y_values)) {
        return *connected;
    }
    return is_feasible_by_max_flow(x_values, y_values, x_offset, y_offset);
}

// Why the traversal gives is_feasible_by_max_flow's answer. Round x and y to
// xr and yr; a used edge has xr_e >= 1. is_feasible_by_max_flow rejects
// exactly when some t != source with y_t > tol has a minimum source–t cut S
// (t in S, source not) with rhs(S)·y_t - x(δ(S)) > tol, where rhs(S) = 1 if
// the path target is in S and 2 otherwise.
//
// - t with yr_t = 0 has y_t <= tol and is skipped there too.
// - Parity: with the degree rows holding for xr and yr, the xr-weighted
//   degrees over S sum to 2·xr(E(S)) + xr(δ(S)) and to an even number plus
//   1 if S holds a path target with yr = 1. So xr(δ(S)) is even unless the
//   path target is in S, and if t is reachable from the source along used
//   edges, xr(δ(S)) >= 1, hence xr(δ(S)) >= rhs(S) for every such S.
// - Reachable t: every used edge has x_e >= xr_e - D (D the largest
//   shortfall), so x(δ(S)) >= xr(δ(S))·(1 - D) >= rhs(S)·(1 - D), while
//   y_t <= 1 + Y (Y the largest excess of a y rounded to 1). The violation is
//   at most rhs(S)·(Y + D) <= 2·(Y + D), which the guard keeps below tol.
// - Unreachable t: S = t's component under used edges is crossed by unused
//   edges only, so the minimum cut is at most M, the unused edges' positive
//   mass, and the violation is at least y_t - M >= 1 - tol - M > tol, with
//   M <= 1/2.
std::optional<bool> SeparationOracle::integral_connectivity(
    std::span<const double> x_values, std::span<const double> y_values) const {
    const auto& graph = prob_.graph();
    const int32_t n = prob_.num_nodes();
    const int32_t source = prob_.source();
    const int32_t target = prob_.target();
    const bool is_tour = prob_.is_tour();

    std::vector<int32_t> xr(static_cast<size_t>(graph.num_edges()), 0);
    double used_shortfall = 0.0;  // D
    double unused_mass = 0.0;     // M
    for (auto e : graph.edges()) {
        const double v = x_values[e];
        const double r = std::round(v);
        // Negated so that a NaN fails the guard too.
        if (!(std::abs(v - r) <= kIncumbentTol) || r < 0.0) {
            return std::nullopt;
        }
        xr[e] = static_cast<int32_t>(r);
        if (xr[e] == 0) {
            unused_mass += std::max(v, 0.0);
        } else {
            used_shortfall = std::max(used_shortfall, r - v);
        }
    }

    std::vector<int8_t> yr(static_cast<size_t>(n), 0);
    double selected_excess = 0.0;  // Y
    for (int32_t i = 0; i < n; ++i) {
        const double v = y_values[i];
        const double r = std::round(v);
        if (!(std::abs(v - r) <= kIncumbentTol) || (r != 0.0 && r != 1.0)) {
            return std::nullopt;
        }
        yr[i] = static_cast<int8_t>(r);
        if (yr[i] == 1) {
            selected_excess = std::max(selected_excess, v - 1.0);
        }
    }

    // Half of tol: room for the rounding in the max-flow's sums.
    if (2.0 * (selected_excess + used_shortfall) > 0.5 * kIncumbentTol || unused_mass > 0.5) {
        return std::nullopt;
    }

    for (int32_t i = 0; i < n; ++i) {
        const bool terminal = !is_tour && (i == source || i == target);
        int64_t degree = 0;
        for (auto e : graph.incident_edges(i)) {
            degree += xr[e];
        }
        if (degree != (terminal ? 1 : 2) * int64_t{yr[i]}) {
            return std::nullopt;
        }
    }

    std::vector<char> reached(static_cast<size_t>(n), 0);
    std::vector<int32_t> stack{source};
    reached[source] = 1;
    while (!stack.empty()) {
        const int32_t u = stack.back();
        stack.pop_back();
        for (auto e : graph.incident_edges(u)) {
            const int32_t w = graph.other_endpoint(e, u);
            if (xr[e] > 0 && reached[w] == 0) {
                reached[w] = 1;
                stack.push_back(w);
            }
        }
    }
    for (int32_t i = 0; i < n; ++i) {
        if (yr[i] == 1 && reached[i] == 0) {
            return false;
        }
    }
    return true;
}

bool SeparationOracle::is_feasible_by_max_flow(std::span<const double> x_values,
                                               std::span<const double> y_values, int32_t x_offset,
                                               int32_t y_offset) const {
    const auto& graph = prob_.graph();
    const int32_t n = prob_.num_nodes();

    digraph_builder builder(n);
    for (auto e : graph.edges()) {
        double xval = x_values[e];
        if (xval > kIncumbentTol) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            builder.add_arc(u, v, xval);
            builder.add_arc(v, u, xval);
        }
    }
    auto [support_graph, capacity] = builder.build();
    gomory_hu_tree flow_tree(support_graph, capacity, prob_.source());

    SeparationContext ctx{
        .problem = prob_,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = x_offset,
        .y_offset = y_offset,
        .tol = kIncumbentTol,
        .flow_tree = &flow_tree,
    };

    SECSeparator sec;
    auto cuts = sec.separate(ctx);
    return cuts.empty();
}

}  // namespace cptp::sep
