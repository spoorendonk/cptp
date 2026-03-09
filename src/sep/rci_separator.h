#pragma once

#include "sep/separator.h"

namespace cptp::sep {

/// Rounded Capacity Inequality separator.
/// For a set S of customers, the RCI is:
///   x(delta(S)) >= 2 * ceil(d(S) / Q)
/// where d(S) = sum of demands in S, Q = vehicle capacity.
class RCISeparator : public Separator {
 public:
  std::string name() const override { return "RCI"; }
  std::vector<Cut> separate(const SeparationContext& ctx) override;
};

}  // namespace cptp::sep
