#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace cptp {

class Problem;

namespace detail {

/// A hyperplane over the y columns, sum_j values[j] * col[indices[j]], branched
/// on through its own integer column z equal to it (#137). `name` labels it.
struct Hyperplane {
    std::vector<int32_t> indices;
    std::vector<double> values;
    std::string name;
};

/// The hyperplanes of `branch_hyper` (pairs, clusters, demand, cardinality or
/// all; "off" or anything else gives none). They are fixed at construction;
/// columns() lists them, each to become one auxiliary integer column.
///
/// The demand-weighted ones (clusters, demand) come per capacity row: row 0
/// under the names it had with one row ("cluster_<j>", "demand"), and each row
/// k >= 1 with 0 < Q_k < kInfiniteCapacity as "cluster_<j>_row<k>" and
/// "demand_row<k>". A row whose demands are not all integer gives none, since
/// z = sum(d_i * y_i) would then not be integer at an integer y; skipped_rows()
/// lists those rows.
class HyperplaneCandidates {
public:
    HyperplaneCandidates(const Problem& problem, int32_t y_offset, std::string mode);

    /// Every hyperplane of the mode, in a fixed order, each once.
    [[nodiscard]] std::vector<Hyperplane> columns() const;

    /// The capacity rows the demand-weighted hyperplanes come from.
    [[nodiscard]] const std::vector<int32_t>& rows() const { return rows_; }

    /// The rows left out of rows() because their demands are not all integer,
    /// when the mode has demand-weighted hyperplanes (clusters, demand, all).
    [[nodiscard]] const std::vector<int32_t>& skipped_rows() const { return skipped_rows_; }

private:
    struct Pair {
        int32_t i, j;
    };

    int32_t y_off_;
    int32_t n_;
    int32_t source_;
    std::string mode_;
    std::vector<int32_t> rows_;
    std::vector<int32_t> skipped_rows_;
    std::vector<std::vector<double>> demands_;  // demands_[r] for rows_[r]
    std::vector<Pair> pairs_;
    std::vector<std::vector<int32_t>> clusters_;

    [[nodiscard]] bool uses(std::string_view m) const { return mode_ == m || mode_ == "all"; }
};

}  // namespace detail
}  // namespace cptp
