#pragma once

#include <condition_variable>
#include <cstdint>
#include <limits>
#include <mutex>
#include <set>

namespace rcspp::model_detail {

// In-process solve admission guard. Finite limits are registered before waiting
// so stricter pending requests constrain later arrivals.
class SolveConcurrencyGuard {
 public:
    explicit SolveConcurrencyGuard(int32_t limit)
        : active_(false),
          requested_limit_(limit > 0 ? limit : kNoCap),
          limited_(requested_limit_ != kNoCap),
          limit_it_(active_limits_.end()) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (limited_) {
            limit_it_ = active_limits_.insert(requested_limit_);
        }
        cv_.wait(lock, [&] {
            return active_solves_ < std::min(requested_limit_, current_cap_locked());
        });
        ++active_solves_;
        active_ = true;
    }

    ~SolveConcurrencyGuard() {
        if (!active_) return;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (limited_) {
                active_limits_.erase(limit_it_);
            }
            --active_solves_;
        }
        cv_.notify_all();
    }

    SolveConcurrencyGuard(const SolveConcurrencyGuard&) = delete;
    SolveConcurrencyGuard& operator=(const SolveConcurrencyGuard&) = delete;

    // Test-only observability helpers.
    static int32_t active_solves_for_testing() {
        std::lock_guard<std::mutex> lock(mutex_);
        return active_solves_;
    }
    static int32_t current_cap_for_testing() {
        std::lock_guard<std::mutex> lock(mutex_);
        return current_cap_locked();
    }

 private:
    static constexpr int32_t kNoCap = std::numeric_limits<int32_t>::max();

    static int32_t current_cap_locked() {
        return active_limits_.empty() ? kNoCap : *active_limits_.begin();
    }

    bool active_;
    int32_t requested_limit_;
    bool limited_;
    std::multiset<int32_t>::iterator limit_it_;

    inline static std::mutex mutex_;
    inline static std::condition_variable cv_;
    inline static int32_t active_solves_ = 0;
    inline static std::multiset<int32_t> active_limits_;
};

}  // namespace rcspp::model_detail

