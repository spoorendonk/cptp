#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Multistar inequality separator.
/// Strengthens capacity cuts using knapsack-like reasoning.
/// With m capacity rows, cuts come from each bounded row's d^k and Q_k as
/// `rows` selects.
class MultistarSeparator : public Separator {
public:
    explicit MultistarSeparator(CapacityRows rows = CapacityRows::kAll) : rows_(rows) {}
    [[nodiscard]] std::string name() const override { return "Multistar"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
};

}  // namespace cptp::sep
