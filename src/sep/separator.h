#pragma once

#include <memory>
#include <string>
#include <vector>

#include "sep/cut.h"
#include "sep/separation_context.h"

namespace cptp::sep {

/// Sentinel capacity value treated as "uncapacitated" by capacity-based
/// separators (GLM, RGLM).  Any capacity >= this threshold disables the cut.
constexpr double kInfiniteCapacity = 1e17;

/// Base class for solver-independent cut separators.
class Separator {
 public:
  virtual ~Separator() = default;
  virtual std::string name() const = 0;
  virtual std::vector<Cut> separate(const SeparationContext& ctx) = 0;
};

}  // namespace cptp::sep
