#include "sep/cut_cleanup.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cptp::sep {

bool drop_near_zero_coefficients(std::vector<int32_t>& idx, std::vector<double>& val, double& rhs,
                                 std::span<const double> lower, std::span<const double> upper,
                                 std::span<const double> lp_point, double min_violation) {
    double max_abs = 0.0;
    for (const double a : val) {
        max_abs = std::max(max_abs, std::abs(a));
    }
    const double threshold = kNearZeroCoefficient * max_abs;

    bool relaxed = false;
    size_t kept = 0;
    for (size_t k = 0; k < idx.size(); ++k) {
        const double a = val[k];
        if (std::abs(a) > threshold) {
            idx[kept] = idx[k];
            val[kept] = a;
            ++kept;
            continue;
        }
        if (a == 0.0) {
            continue;
        }
        const auto col = static_cast<size_t>(idx[k]);
        const double bound = a > 0.0 ? lower[col] : upper[col];
        if (!std::isfinite(bound)) {
            return false;
        }
        rhs -= a * bound;
        relaxed = true;
    }
    idx.resize(kept);
    val.resize(kept);
    if (kept == 0) {
        return false;
    }
    if (!relaxed) {
        return true;
    }
    double activity = 0.0;
    for (size_t k = 0; k < kept; ++k) {
        activity += val[k] * lp_point[static_cast<size_t>(idx[k])];
    }
    return activity - rhs > min_violation;
}

}  // namespace cptp::sep
