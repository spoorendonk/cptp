#pragma once

#include <cstdint>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

namespace rcspp::model {

struct AsyncIncumbentSnapshot {
    std::vector<double> col_values;
    double objective = std::numeric_limits<double>::infinity();
    uint64_t version = 0;
};

class AsyncIncumbentStore {
 public:
    AsyncIncumbentSnapshot snapshot() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return snapshot_;
    }

    bool publish_if_better(std::vector<double> col_values, double objective) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (col_values.empty()) return false;
        if (!(objective + 1e-9 < snapshot_.objective)) return false;
        snapshot_.objective = objective;
        snapshot_.col_values = std::move(col_values);
        snapshot_.version++;
        return true;
    }

 private:
    mutable std::mutex mutex_;
    AsyncIncumbentSnapshot snapshot_;
};

}  // namespace rcspp::model
