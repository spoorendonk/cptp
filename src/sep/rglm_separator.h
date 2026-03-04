#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded GLM (RGLM) inequality separator (eq 24, Jepsen et al. 2014).
/// Strengthens GLM/multistar cuts using ceiling-based rounding.
class RGLMSeparator : public Separator {
 public:
    std::string name() const override { return "RGLM"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace cptp::sep
