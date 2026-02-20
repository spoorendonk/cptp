#pragma once

#include <chrono>

namespace cptp {

class Timer {
 public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    void reset() { start_ = std::chrono::steady_clock::now(); }

    double elapsed_seconds() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - start_).count();
    }

 private:
    std::chrono::steady_clock::time_point start_;
};

}  // namespace cptp
