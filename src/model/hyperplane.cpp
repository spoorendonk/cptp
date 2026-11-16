#include "model/hyperplane.h"

#include "core/problem.h"
#include "sep/separator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace cptp::detail {

HyperplaneCandidates::HyperplaneCandidates(const Problem& problem, int32_t y_offset,
                                           std::string mode)
    : y_off_(y_offset), n_(problem.num_nodes()), source_(problem.source()), mode_(std::move(mode)) {
    const int32_t n = n_;
    const auto& graph = problem.graph();

    // Integer demands give an integer z at an integer y (same 1e-9 test as
    // the bridge's integer-objective check).
    auto integer_demands = [&](int32_t k) {
        return std::ranges::all_of(problem.demands(k),
                                   [](double d) { return std::abs(d - std::round(d)) < 1e-9; });
    };
    for (int32_t k = 0; k < problem.num_capacities(); ++k) {
        const double q = problem.capacity(k);
        if (k > 0 && (q <= 0 || q >= sep::kInfiniteCapacity)) {
            continue;
        }
        if (integer_demands(k)) {
            rows_.push_back(k);
        } else {
            skipped_rows_.push_back(k);
        }
    }
    for (int32_t k : rows_) {
        demands_.push_back(problem.demands(k));
    }
    const auto& demands = problem.demands(0);

    // Ryan-Foster pairs: for each non-depot node, pair with nearest neighbor.
    if (uses("pairs")) {
        std::vector<bool> paired(n, false);
        for (int32_t i = 0; i < n; ++i) {
            if (i == source_) {
                continue;
            }
            double best_cost = std::numeric_limits<double>::infinity();
            int32_t best_j = -1;
            for (auto e : graph.incident_edges(i)) {
                int32_t u = graph.edge_source(e);
                int32_t v = graph.edge_target(e);
                int32_t nb = (u == i) ? v : u;
                if (nb == source_) {
                    continue;
                }
                double c = problem.edge_cost(e);
                if (c < best_cost) {
                    best_cost = c;
                    best_j = nb;
                }
            }
            if (best_j >= 0) {
                // Deduplicate: only add (min, max) pairs.
                int32_t lo = std::min(i, best_j);
                int32_t hi = std::max(i, best_j);
                if (!paired[lo] && !paired[hi]) {
                    pairs_.push_back({lo, hi});
                    paired[lo] = true;
                    paired[hi] = true;
                }
            }
        }
    }

    // Greedy nearest-neighbor clusters of size ~4, seeded by row 0's demand.
    // The clusters are geometric; each row weighs them with its own demands.
    if (uses("clusters")) {
        constexpr int32_t kClusterSize = 4;
        std::vector<bool> assigned(n, false);
        assigned[source_] = true;

        // Sort nodes by demand (descending) to seed clusters.
        std::vector<int32_t> nodes_by_demand;
        for (int32_t i = 0; i < n; ++i) {
            if (i != source_) {
                nodes_by_demand.push_back(i);
            }
        }
        std::ranges::sort(nodes_by_demand,
                          [&](int32_t a, int32_t b) { return demands[a] > demands[b]; });

        for (int32_t seed : nodes_by_demand) {
            if (assigned[seed]) {
                continue;
            }
            std::vector<int32_t> cluster = {seed};
            assigned[seed] = true;

            // Add nearest unassigned neighbors.
            while (static_cast<int32_t>(cluster.size()) < kClusterSize) {
                double best_cost = std::numeric_limits<double>::infinity();
                int32_t best_nb = -1;
                for (int32_t ci : cluster) {
                    for (auto e : graph.incident_edges(ci)) {
                        int32_t u = graph.edge_source(e);
                        int32_t v = graph.edge_target(e);
                        int32_t nb = (u == ci) ? v : u;
                        if (assigned[nb]) {
                            continue;
                        }
                        double c = problem.edge_cost(e);
                        if (c < best_cost) {
                            best_cost = c;
                            best_nb = nb;
                        }
                    }
                }
                if (best_nb < 0) {
                    break;
                }
                cluster.push_back(best_nb);
                assigned[best_nb] = true;
            }
            if (cluster.size() >= 2) {
                clusters_.push_back(std::move(cluster));
            }
        }
    }
}

std::vector<Hyperplane> HyperplaneCandidates::columns() const {
    std::vector<Hyperplane> result;

    // Row 0 keeps its one-row names; row k >= 1 gets "_row<k>".
    auto suffix = [this](size_t r) {
        return rows_[r] == 0 ? std::string{} : "_row" + std::to_string(rows_[r]);
    };

    if (uses("pairs")) {
        for (const auto& [pi, pj] : pairs_) {
            Hyperplane c;
            c.indices = {y_off_ + pi, y_off_ + pj};
            c.values = {1.0, 1.0};
            c.name = "pair_" + std::to_string(pi) + "_" + std::to_string(pj);
            result.push_back(std::move(c));
        }
    }
    if (uses("clusters")) {
        for (size_t r = 0; r < rows_.size(); ++r) {
            const auto& demands = demands_[r];
            for (size_t k = 0; k < clusters_.size(); ++k) {
                // Cluster demand hyperplane: sum(demand_i * y_i) for i in cluster.
                Hyperplane c;
                c.name = "cluster_" + std::to_string(k) + suffix(r);
                for (int32_t i : clusters_[k]) {
                    if (demands[i] > 0) {
                        c.indices.push_back(y_off_ + i);
                        c.values.push_back(demands[i]);
                    }
                }
                if (!c.indices.empty()) {
                    result.push_back(std::move(c));
                }
            }
        }
    }
    if (uses("demand")) {
        for (size_t r = 0; r < rows_.size(); ++r) {
            const auto& demands = demands_[r];
            Hyperplane c;
            c.name = "demand" + suffix(r);
            for (int32_t i = 0; i < n_; ++i) {
                if (i == source_) {
                    continue;  // skip depot
                }
                if (demands[i] > 0) {
                    c.indices.push_back(y_off_ + i);
                    c.values.push_back(demands[i]);
                }
            }
            if (!c.indices.empty()) {
                result.push_back(std::move(c));
            }
        }
    }
    if (uses("cardinality")) {
        Hyperplane c;
        c.name = "cardinality";
        for (int32_t i = 0; i < n_; ++i) {
            if (i == source_) {
                continue;  // skip depot
            }
            c.indices.push_back(y_off_ + i);
            c.values.push_back(1.0);
        }
        if (!c.indices.empty()) {
            result.push_back(std::move(c));
        }
    }
    return result;
}

}  // namespace cptp::detail
