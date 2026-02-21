#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace rcspp {

/// A resource constraint on the problem (e.g., capacity).
struct Resource {
    enum class Type { Capacity };

    Type type;
    std::vector<double> node_consumption;  // per-node resource usage
    double limit;                          // upper bound
};

}  // namespace rcspp
