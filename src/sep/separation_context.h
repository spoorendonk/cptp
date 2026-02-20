#pragma once

#include <cstdint>
#include <span>

namespace cptp {
class Problem;
}

namespace cptp::sep {

/// Solver-independent context passed to each separator.
struct SeparationContext {
    const Problem& problem;
    std::span<const double> x_values;  // edge variables
    std::span<const double> y_values;  // node variables
    int32_t x_offset;  // LP column offset for edge vars
    int32_t y_offset;  // LP column offset for node vars
    double tol = 1e-6;
};

}  // namespace cptp::sep
