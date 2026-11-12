#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cptp {

class Problem;

namespace detail {

/// A hyperplane branching candidate over the y columns: sum(values * x[indices]).
/// `name` keys its pseudocosts, so the same hyperplane keeps the same name.
struct Hyperplane {
    std::vector<int32_t> indices;
    std::vector<double> values;
    std::string name;
};

/// The candidates of `branch_hyper` (pairs, clusters, demand, cardinality or
/// all; "off" or anything else gives none). The pairs and clusters are fixed
/// at construction; operator() lists the candidates at an LP solution.
///
/// The demand-weighted ones (clusters, demand) come per capacity row: row 0
/// always, under the names it had with one row ("cluster_<j>", "demand"), and
/// each row k >= 1 with 0 < Q_k < kInfiniteCapacity as "cluster_<j>_row<k>" and
/// "demand_row<k>", so each row's pseudocosts are kept apart. Every such row,
/// not only the most fractional: which one to branch on is the pseudocost and
/// strong-branching selection's job, and a candidate that comes and goes with
/// the LP would keep its pseudocosts from building up (#100).
class HyperplaneCandidates {
public:
    HyperplaneCandidates(const Problem& problem, int32_t y_offset, std::string mode);

    [[nodiscard]] std::vector<Hyperplane> operator()(std::span<const double> col_value) const;

    /// The capacity rows the demand-weighted candidates come from.
    [[nodiscard]] const std::vector<int32_t>& rows() const { return rows_; }

private:
    struct Pair {
        int32_t i, j;
    };

    int32_t y_off_;
    int32_t n_;
    int32_t source_;
    std::string mode_;
    std::vector<int32_t> rows_;
    std::vector<std::vector<double>> demands_;  // demands_[r] for rows_[r]
    std::vector<Pair> pairs_;
    std::vector<std::vector<int32_t>> clusters_;

    [[nodiscard]] bool uses(std::string_view m) const { return mode_ == m || mode_ == "all"; }
};

}  // namespace detail
}  // namespace cptp
