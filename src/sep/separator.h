#pragma once

#include <memory>
#include <string>
#include <vector>

#include "sep/cut.h"
#include "sep/separation_context.h"

namespace cptp::sep {

/// Base class for solver-independent cut separators.
class Separator {
 public:
  virtual ~Separator() = default;
  virtual std::string name() const = 0;
  virtual std::vector<Cut> separate(const SeparationContext& ctx) = 0;
};

}  // namespace cptp::sep
