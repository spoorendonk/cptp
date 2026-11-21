#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Multistar inequality separator.
/// Strengthens capacity cuts using knapsack-like reasoning.
/// With m capacity rows, cuts come from each bounded row's d^k and Q_k as
/// `rows` selects, Q_k less the terminals' demand and the terminals weighted
/// 0 (residual_capacity, #170).
///
/// Each target node's set S is the smallest depot-target min-cut side the
/// Gomory-Hu tree gives. With `largest_side` it also tries the largest one
/// (#168): one depot-to-target max-flow on the support graph per target, S
/// being the nodes the residual does not reach from the depot that are
/// connected to the target within that set (so unvisited customers and
/// components cut off from both ends stay out). Each S is cut once. It needs
/// the context's support graph and is skipped without it.
class MultistarSeparator : public Separator {
public:
    explicit MultistarSeparator(CapacityRows rows = CapacityRows::kAll, bool largest_side = false)
        : rows_(rows), largest_side_(largest_side) {}
    [[nodiscard]] std::string name() const override { return "Multistar"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
    bool largest_side_;
};

}  // namespace cptp::sep
