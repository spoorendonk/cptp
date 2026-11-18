#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>

namespace cptp::parallel {

namespace detail {
/// Threads parallel_for has spawned that are still running, and the most that
/// ran at once since the last reset_peak_spawned(): lets a test check that a
/// solve kept to its thread budget.
inline std::atomic<int> live_spawned{0};
inline std::atomic<int> peak_spawned{0};

/// Counts one spawned thread for its lifetime.
class SpawnedThreadScope {
public:
    SpawnedThreadScope() {
        const int live = live_spawned.fetch_add(1, std::memory_order_relaxed) + 1;
        int peak = peak_spawned.load(std::memory_order_relaxed);
        while (live > peak &&
               !peak_spawned.compare_exchange_weak(peak, live, std::memory_order_relaxed)) {
        }
    }
    ~SpawnedThreadScope() { live_spawned.fetch_sub(1, std::memory_order_relaxed); }
    SpawnedThreadScope(const SpawnedThreadScope&) = delete;
    SpawnedThreadScope& operator=(const SpawnedThreadScope&) = delete;
    SpawnedThreadScope(SpawnedThreadScope&&) = delete;
    SpawnedThreadScope& operator=(SpawnedThreadScope&&) = delete;
};
}  // namespace detail

/// Most threads parallel_for had spawned at once since the last reset.
inline int peak_spawned_threads() {
    return detail::peak_spawned.load(std::memory_order_relaxed);
}

/// Restart the peak count from the threads running now.
inline void reset_peak_spawned_threads() {
    detail::peak_spawned.store(detail::live_spawned.load(std::memory_order_relaxed),
                               std::memory_order_relaxed);
}

/// Split [begin, end) across at most num_threads threads, the calling thread
/// among them: it spawns num_threads - 1 and runs the last chunk itself.
/// num_threads <= 1 runs everything on the calling thread.
template <typename F>
void parallel_for(int begin, int end, F&& f, int num_threads) {
    if (begin >= end) {
        return;
    }
    const int range = end - begin;
    num_threads = std::clamp(num_threads, 1, range);
    if (num_threads == 1) {
        for (int i = begin; i < end; ++i) {
            f(i);
        }
        return;
    }

    std::vector<std::jthread> threads;
    threads.reserve(static_cast<size_t>(num_threads - 1));
    const int chunk = range / num_threads;
    const int remainder = range % num_threads;

    int lo = begin;
    for (int t = 0; t < num_threads - 1; ++t) {
        const int hi = lo + chunk + (t < remainder ? 1 : 0);
        threads.emplace_back([&f, lo, hi] {
            const detail::SpawnedThreadScope counted;
            for (int i = lo; i < hi; ++i) {
                f(i);
            }
        });
        lo = hi;
    }
    for (int i = lo; i < end; ++i) {
        f(i);
    }
    // jthreads auto-join on destruction
}

/// Spawn tasks that auto-join on destruction (thin jthread wrapper). A serial
/// group runs each task on the calling thread inside run().
class task_group {
    std::vector<std::jthread> threads_;
    bool serial_ = false;

public:
    task_group() = default;
    explicit task_group(bool serial) : serial_(serial) {}

    template <typename F>
    void run(F&& f) {
        if (serial_) {
            std::forward<F>(f)();
            return;
        }
        threads_.emplace_back([fn = std::forward<F>(f)] { fn(); });
    }

    void wait() {
        for (auto& t : threads_) {
            if (t.joinable()) {
                t.join();
            }
        }
        threads_.clear();
    }

    ~task_group() { wait(); }
};

}  // namespace cptp::parallel
