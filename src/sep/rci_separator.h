#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded Capacity Inequality separator.
/// For a set S of customers, the RCI is:
///   x(delta(S)) >= 2 * ceil(d(S) / Q)
/// where d(S) = sum of demands in S, Q = vehicle capacity.
/// With m capacity rows, cuts are built from each row's d^k and Q_k as
/// `rows` selects (unbounded rows, Q_k >= kInfiniteCapacity, give none).
/// With `surrogate`, one more row next to them: the surrogate row of all
/// bounded rows (surrogate_row, #112). Under kMostViolated it competes with
/// the others; under kAll and kFirst its cuts are added to theirs.
class RCISeparator : public Separator {
public:
    explicit RCISeparator(CapacityRows rows = CapacityRows::kAll, bool surrogate = false)
        : rows_(rows), surrogate_(surrogate) {}
    [[nodiscard]] std::string name() const override { return "RCI"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
    bool surrogate_;
};

}  // namespace cptp::sep
