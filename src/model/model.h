#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/problem.h"
#include "core/solution.h"
#include "model/highs_bridge.h"
#include "model/resource.h"

namespace cptp {

/// User-facing model: build a CPTP instance and solve it.
class Model {
 public:
    Model();

    /// Set problem directly (e.g., from io::load).
    void set_problem(Problem prob);

    void set_graph(int32_t num_nodes,
                   std::span<const Edge> edges,
                   std::span<const double> edge_costs);

    void set_source(int32_t source);
    void set_target(int32_t target);
    void set_depot(int32_t depot);  // sets source = target = depot (tour)
    void set_profits(std::span<const double> profits);

    void add_capacity_resource(std::span<const double> demands, double limit);

    SolveResult solve(const SolverOptions& options = {});

    const Problem& problem() const { return problem_; }

 private:
    void build_problem();

    int32_t num_nodes_ = 0;
    int32_t source_ = 0;
    int32_t target_ = 0;
    std::vector<Edge> edges_;
    std::vector<double> edge_costs_;
    std::vector<double> profits_;
    std::vector<Resource> resources_;
    Problem problem_;
    bool built_ = false;
};

}  // namespace cptp
