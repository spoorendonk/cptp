#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded GLM (RGLM) inequality separator (eq 24, Jepsen et al. 2014).
/// Strengthens GLM/multistar cuts using ceiling-based rounding.
/// With m capacity rows, cuts come from each bounded row's d^k and Q_k as
/// `rows` selects.
/// With `surrogate`, one more row next to them: the surrogate row of all
/// bounded rows (surrogate_row, #112). Under kMostViolated it competes with
/// the others; under kAll and kFirst its cuts are added to theirs.
class RGLMSeparator : public Separator {
public:
    explicit RGLMSeparator(CapacityRows rows = CapacityRows::kAll, bool surrogate = false)
        : rows_(rows), surrogate_(surrogate) {}
    [[nodiscard]] std::string name() const override { return "RGLM"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
    bool surrogate_;
};

}  // namespace cptp::sep
