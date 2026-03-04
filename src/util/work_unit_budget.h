#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>

namespace cptp {

/// Thread-safe work-unit budget.
/// limit <= 0 means "count only, no cap".
class WorkUnitBudget {
 public:
    explicit WorkUnitBudget(int64_t limit = 0)
        : limit_(std::max<int64_t>(0, limit)) {}

    bool capped() const { return limit_ > 0; }
    int64_t limit() const { return limit_; }

    int64_t used() const {
        return used_.load(std::memory_order_relaxed);
    }

    int64_t remaining() const {
        if (!capped()) return -1;
        const int64_t rem = limit_ - used();
        return std::max<int64_t>(0, rem);
    }

    /// Try to consume a fixed number of work units.
    /// Returns false only when capped and insufficient units remain.
    bool try_consume(int64_t units = 1) {
        if (units <= 0) return true;
        if (!capped()) {
            used_.fetch_add(units, std::memory_order_relaxed);
            return true;
        }

        int64_t cur = used_.load(std::memory_order_relaxed);
        while (true) {
            if (cur >= limit_) return false;
            const int64_t grant = std::min<int64_t>(units, limit_ - cur);
            if (grant < units) return false;
            if (used_.compare_exchange_weak(
                    cur, cur + units,
                    std::memory_order_relaxed,
                    std::memory_order_relaxed)) {
                return true;
            }
        }
    }

    /// Reserve up to `units` in a single atomic step.
    /// Returns the granted amount in [0, units].
    int64_t reserve_up_to(int64_t units) {
        if (units <= 0) return 0;
        if (!capped()) {
            used_.fetch_add(units, std::memory_order_relaxed);
            return units;
        }

        int64_t cur = used_.load(std::memory_order_relaxed);
        while (true) {
            if (cur >= limit_) return 0;
            const int64_t grant = std::min<int64_t>(units, limit_ - cur);
            if (used_.compare_exchange_weak(
                    cur, cur + grant,
                    std::memory_order_relaxed,
                    std::memory_order_relaxed)) {
                return grant;
            }
        }
    }

 private:
    int64_t limit_ = 0;
    std::atomic<int64_t> used_{0};
};

}  // namespace cptp
