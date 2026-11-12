#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded GLM (RGLM) inequality separator (eq 24, Jepsen et al. 2014).
/// Strengthens GLM/multistar cuts using ceiling-based rounding.
/// With m capacity rows, cuts come from each bounded row's d^k and Q_k as
/// `rows` selects.
class RGLMSeparator : public Separator {
public:
    explicit RGLMSeparator(CapacityRows rows = CapacityRows::kAll) : rows_(rows) {}
    [[nodiscard]] std::string name() const override { return "RGLM"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;

private:
    CapacityRows rows_;
};

}  // namespace cptp::sep
