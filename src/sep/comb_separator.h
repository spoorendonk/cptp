#pragma once

#include "sep/separator.h"

namespace rcspp::sep {

/// Comb inequality separator.
class CombSeparator : public Separator {
 public:
    std::string name() const override { return "Comb"; }
    std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace rcspp::sep
