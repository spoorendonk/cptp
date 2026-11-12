#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded Capacity Inequality separator.
/// For a set S of customers, the RCI is:
///   x(delta(S)) >= 2 * ceil(d(S) / Q)
/// where d(S) = sum of demands in S, Q = vehicle capacity.
/// With m capacity rows, cuts are built from each row's d^k and Q_k as
/// `rows` selects (unbounded rows, Q_k >= kInfiniteCapacity, give none).
class RCISeparator : public Separator {
public:
    explicit RCISeparator(CapacityRows rows = CapacityRows::kAll) : rows_(rows) {}
    [[nodiscard]] std::string name() const override { return "RCI"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
};

}  // namespace cptp::sep
