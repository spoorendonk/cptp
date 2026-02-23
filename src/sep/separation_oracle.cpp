#include "sep/separation_oracle.h"

#include <algorithm>
#include <chrono>
#include <vector>

#include <tbb/task_group.h>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "sep/comb_separator.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/sec_separator.h"

namespace rcspp::sep {

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

std::vector<Cut> SeparationOracle::separate(
    std::span<const double> x_values,
    std::span<const double> y_values,
    int32_t x_offset,
    int32_t y_offset,
    double tol) const {

    const auto& graph = prob_.graph();
    const int32_t n = prob_.num_nodes();
    const int32_t m = prob_.num_edges();

    // Build support graph from fractional LP values
    constexpr double graph_tol = 1e-6;
    digraph_builder builder(n);
    for (auto e : graph.edges()) {
        double xval = x_values[e];
        if (xval > graph_tol) {
            int32_t u = graph.edge_source(e);
            int32_t v = graph.edge_target(e);
            builder.add_arc(u, v, xval);
            builder.add_arc(v, u, xval);
        }
    }
    auto [support_graph, capacity] = builder.build();

    // Compute Gomory-Hu tree (n-1 max-flows via Gusfield)
    gomory_hu_tree flow_tree(support_graph, capacity, prob_.source());

    SeparationContext ctx{
        .problem = prob_,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = x_offset,
        .y_offset = y_offset,
        .tol = tol,
        .flow_tree = &flow_tree,
    };

    // Run all separators in parallel
    std::vector<std::vector<Cut>> results(separators_.size());
    tbb::task_group tg;
    for (size_t i = 0; i < separators_.size(); ++i) {
        tg.run([&, i] {
            results[i] = separators_[i]->separate(ctx);
        });
    }
    tg.wait();

    // Collect, sort by violation, limit per separator
    std::vector<Cut> all_cuts;
    for (size_t i = 0; i < separators_.size(); ++i) {
        auto& cuts = results[i];
        std::sort(cuts.begin(), cuts.end(),
            [](const Cut& a, const Cut& b) {
                return a.violation > b.violation;
            });
        int32_t limit = (max_cuts_per_sep_ > 0)
            ? std::min(static_cast<int32_t>(cuts.size()), max_cuts_per_sep_)
            : static_cast<int32_t>(cuts.size());
        for (int32_t j = 0; j < limit; ++j) {
            all_cuts.push_back(std::move(cuts[j]));
        }
    }

    // Sort combined results by violation
    std::sort(all_cuts.begin(), all_cuts.end(),
        [](const Cut& a, const Cut& b) {
            return a.violation > b.violation;
        });

    return all_cuts;
}

bool SeparationOracle::is_feasible(
    std::span<const double> x_values,
    std::span<const double> y_values,
    int32_t x_offset,
    int32_t y_offset) const {

    const auto& graph = prob_.graph();
    const int32_t n = prob_.num_nodes();

    constexpr double graph_tol = 1e-6;
    digraph_builder builder(n);
    for (auto e : graph.edges()) {
        double xval = x_values[e];
        if (xval > graph_tol) {
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
        .tol = 1e-6,
        .flow_tree = &flow_tree,
    };

    SECSeparator sec;
    auto cuts = sec.separate(ctx);
    return cuts.empty();
}

}  // namespace rcspp::sep
