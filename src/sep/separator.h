#pragma once

#include "sep/cut.h"
#include "sep/separation_context.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cptp::sep {

/// Sentinel capacity value treated as "uncapacitated" by capacity-based
/// separators (GLM, RGLM).  Any capacity >= this threshold disables the cut.
constexpr double kInfiniteCapacity = 1e17;

/// Which capacity rows a capacity-based separator builds cuts from (#91).
/// Each row k gives cuts from d^k and Q_k that are valid on their own.
enum class CapacityRows {
    kAll,           ///< every bounded row, each round
    kMostViolated,  ///< per candidate set, only the row whose cut is most violated
    kFirst,         ///< row 0 only (ablation; the single-row behaviour)
};

/// Parse the `capacity_cut_rows` option: "all", "most_violated" or "first".
inline CapacityRows parse_capacity_rows(std::string_view value) {
    if (value == "all") {
        return CapacityRows::kAll;
    }
    if (value == "most_violated") {
        return CapacityRows::kMostViolated;
    }
    if (value == "first") {
        return CapacityRows::kFirst;
    }
    throw std::invalid_argument("capacity_cut_rows must be all, most_violated or first, not " +
                                std::string(value));
}

/// Base class for solver-independent cut separators.
class Separator {
public:
    virtual ~Separator() = default;
    [[nodiscard]] virtual std::string name() const = 0;
    virtual std::vector<Cut> separate(const SeparationContext& ctx) = 0;
};

}  // namespace cptp::sep
