#pragma once

#include "core/problem.h"
#include "sep/separator.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cptp::sep {

/// One surrogate capacity row: demand s(v) = sum_k d_k(v) / Q_k against
/// capacity S = the number of rows summed, over the rows with
/// 0 < Q_k < kInfiniteCapacity. A path or tour within every row is within it
/// (sum of terms <= 1 each), so it relaxes the rows: labeling it gives a bound
/// (#102), and a capacity cut from it is valid (#112). A row left out only
/// loosens it further; Q_k = 0 is left out because 1/Q_k is undefined. S
/// carries a relative slack of 1e-9, so rounding in the scaling cannot cut
/// off a path that fits. Empty demand when fewer than two rows qualify: then
/// the surrogate is no stronger than one row.
struct SurrogateRow {
    std::vector<double> demand;
    double capacity = 0.0;
};

inline SurrogateRow surrogate_row(const Problem& prob) {
    SurrogateRow row;
    int32_t used = 0;
    for (int32_t k = 0; k < prob.num_capacities(); ++k) {
        const double q = prob.capacity(k);
        if (q <= 0.0 || q >= kInfiniteCapacity) {
            continue;
        }
        if (row.demand.empty()) {
            row.demand.assign(static_cast<size_t>(prob.num_nodes()), 0.0);
        }
        for (int32_t v = 0; v < prob.num_nodes(); ++v) {
            row.demand[static_cast<size_t>(v)] += prob.demand(k, v) / q;
        }
        ++used;
    }
    if (used < 2) {
        return {};
    }
    row.capacity = static_cast<double>(used) * (1.0 + 1e-9);
    return row;
}

/// A capacity row a capacity cut is built from: demand per node and capacity.
struct CapacityRowView {
    std::span<const double> demand;
    double capacity;
};

/// The rows a capacity separator builds cuts from: the rows `rows` selects
/// with 0 < Q_k < kInfiniteCapacity, then `surrogate` when its demand is not
/// empty. The views point into `prob` and `surrogate`, which must outlive
/// them.
inline std::vector<CapacityRowView> capacity_row_views(const Problem& prob, CapacityRows rows,
                                                       const SurrogateRow& surrogate) {
    std::vector<CapacityRowView> views;
    const int32_t num_rows = rows == CapacityRows::kFirst ? 1 : prob.num_capacities();
    for (int32_t k = 0; k < num_rows; ++k) {
        const double q = prob.capacity(k);
        if (q > 0 && q < kInfiniteCapacity) {
            views.push_back({prob.demands(k), q});
        }
    }
    if (!surrogate.demand.empty()) {
        views.push_back({surrogate.demand, surrogate.capacity});
    }
    return views;
}

}  // namespace cptp::sep
