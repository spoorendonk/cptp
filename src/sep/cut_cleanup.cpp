#include "sep/cut_cleanup.h"

#include <cmath>
#include <cstddef>

namespace cptp::sep {

std::optional<double> drop_near_zero_coefficients(std::vector<int32_t>& idx,
                                                  std::vector<double>& val, double& rhs,
                                                  std::span<const double> lower,
                                                  std::span<const double> upper,
                                                  std::span<const double> lp_point) {
    double lost = 0.0;
    size_t kept = 0;
    for (size_t k = 0; k < idx.size(); ++k) {
        const double a = val[k];
        if (!std::isfinite(a)) {
            return std::nullopt;
        }
        if (std::abs(a) > kNearZeroCoefficient) {
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
            return std::nullopt;
        }
        rhs -= a * bound;
        lost += a * (lp_point[col] - bound);
    }
    idx.resize(kept);
    val.resize(kept);
    if (kept == 0) {
        return std::nullopt;
    }
    return lost;
}

}  // namespace cptp::sep
