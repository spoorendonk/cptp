#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>

namespace rcspp::preprocess {

struct BoundSnapshot {
    std::vector<double> fwd;
    std::vector<double> bwd;
    double correction = 0.0;
    int32_t ng_size = 1;
    uint64_t version = 0;
};

/// Thread-safe, versioned bounds storage for asynchronous DSSR updates.
class SharedBoundsStore {
 public:
    explicit SharedBoundsStore(BoundSnapshot initial)
        : snapshot_(std::move(initial)) {}

    BoundSnapshot snapshot() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return snapshot_;
    }

    /// Publish tighter bounds (element-wise max on finite values).
    /// Returns true if snapshot improved.
    bool publish_tightening(const std::vector<double>& fwd,
                            const std::vector<double>& bwd,
                            int32_t ng_size) {
        std::lock_guard<std::mutex> lock(mutex_);
        bool changed = false;
        changed |= tighten_vector(snapshot_.fwd, fwd);
        changed |= tighten_vector(snapshot_.bwd, bwd);
        if (ng_size > snapshot_.ng_size) {
            snapshot_.ng_size = ng_size;
            changed = true;
        }
        if (changed) snapshot_.version++;
        return changed;
    }

 private:
    static bool tighten_vector(std::vector<double>& dst,
                               const std::vector<double>& src) {
        if (dst.size() != src.size()) return false;
        constexpr double neg_inf = -std::numeric_limits<double>::infinity();
        bool changed = false;
        for (size_t i = 0; i < dst.size(); ++i) {
            const double cur = dst[i];
            const double nxt = src[i];
            if (!std::isfinite(nxt)) continue;
            if (!std::isfinite(cur)) {
                dst[i] = nxt;
                changed = true;
                continue;
            }
            // Lower bounds should tighten upward.
            const double tightened = std::max(cur, nxt);
            if (tightened > cur && tightened > neg_inf) {
                dst[i] = tightened;
                changed = true;
            }
        }
        return changed;
    }

    mutable std::mutex mutex_;
    BoundSnapshot snapshot_;
};

}  // namespace rcspp::preprocess
