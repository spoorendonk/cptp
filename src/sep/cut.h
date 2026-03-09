#pragma once

#include <cstdint>
#include <vector>

namespace cptp::sep {

/// Sparse representation of a single linear cut: a^T x <= rhs.
struct Cut {
  std::vector<int32_t> indices;
  std::vector<double> values;
  double rhs;
  double violation = 0.0;  // how much the current LP violates this cut

  int32_t size() const { return static_cast<int32_t>(indices.size()); }
};

}  // namespace cptp::sep
