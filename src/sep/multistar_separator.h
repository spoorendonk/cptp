#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Multistar inequality separator.
/// Strengthens capacity cuts using knapsack-like reasoning.
/// With m capacity rows, cuts come from each bounded row's d^k and Q_k as
/// `rows` selects.
///
/// Each target node's set S is the smallest depot-target min-cut side the
/// Gomory-Hu tree gives. With `largest_side` it also tries the largest one
/// (#168): one depot-to-target max-flow on the support graph per target, S
/// being the nodes the residual does not reach from the depot. Among tied
/// min cuts that side gives the most violated GLM. It needs the context's
/// support graph and is skipped without it.
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
