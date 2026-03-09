#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Multistar inequality separator.
/// Strengthens capacity cuts using knapsack-like reasoning.
class MultistarSeparator : public Separator {
 public:
  std::string name() const override { return "Multistar"; }
  std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace cptp::sep
