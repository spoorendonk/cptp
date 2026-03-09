#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Subtour Elimination Constraint separator.
/// Uses Dinitz max-flow on the fractional support graph to find
/// violated min-cuts between the depot and each customer node.
class SECSeparator : public Separator {
 public:
  std::string name() const override { return "SEC"; }
  std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace cptp::sep
