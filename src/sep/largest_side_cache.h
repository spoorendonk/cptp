#pragma once

#include "core/dinitz.h"

#include <atomic>
#include <cstdint>
#include <mutex>
#include <vector>

namespace cptp::sep {

/// Each target's largest depot–target min-cut sink side on one support graph
/// (largest_min_cut_sink_side, connected to the target as #168 has it), run on
/// a target's first request and kept for the round. The multistar
/// separator's largest-side option and the cross-row covers ask for the same
/// targets, and a round runs its separators in parallel, so get() may be
/// called from several threads at once: each target's max-flow runs once,
/// and a second caller waits for it (#175).
class LargestSideCache {
public:
    /// `g` and `cap` must outlive the cache.
    LargestSideCache(const digraph& g, const std::vector<double>& cap, int32_t source)
        : g_(g),
          cap_(cap),
          source_(source),
          once_(static_cast<size_t>(g.num_nodes())),
          sides_(static_cast<size_t>(g.num_nodes())) {}

    /// The largest sink side of a minimum `source`–`target` cut.
    const LargestSinkSide& get(int32_t target) {
        std::call_once(once_[target], [&] {
            sides_[target] = largest_min_cut_sink_side(g_, cap_, source_, target);
            max_flows_.fetch_add(1, std::memory_order_relaxed);
        });
        return sides_[target];
    }

    /// Max-flows run so far, one per target asked for.
    [[nodiscard]] int64_t max_flows() const { return max_flows_.load(std::memory_order_relaxed); }

private:
    const digraph& g_;
    const std::vector<double>& cap_;
    int32_t source_;
    std::vector<std::once_flag> once_;
    std::vector<LargestSinkSide> sides_;
    std::atomic<int64_t> max_flows_{0};
};

}  // namespace cptp::sep
