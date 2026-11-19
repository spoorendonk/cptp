#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>

namespace cptp::parallel {

namespace detail {
/// Threads parallel_for has spawned that are still running, with the outside
/// threads count_busy_threads reports, and the most that ran at once since
/// the last reset_peak_spawned(): lets a test check that a solve kept to its
/// thread budget.
inline std::atomic<int> live_spawned{0};
inline std::atomic<int> peak_spawned{0};

inline void add_live(int delta) {
    const int live = live_spawned.fetch_add(delta, std::memory_order_relaxed) + delta;
    int peak = peak_spawned.load(std::memory_order_relaxed);
    while (live > peak &&
           !peak_spawned.compare_exchange_weak(peak, live, std::memory_order_relaxed)) {
    }
}

/// Counts one spawned thread for its lifetime.
class SpawnedThreadScope {
public:
    SpawnedThreadScope() { add_live(1); }
    ~SpawnedThreadScope() { add_live(-1); }
    SpawnedThreadScope(const SpawnedThreadScope&) = delete;
    SpawnedThreadScope& operator=(const SpawnedThreadScope&) = delete;
    SpawnedThreadScope(SpawnedThreadScope&&) = delete;
    SpawnedThreadScope& operator=(SpawnedThreadScope&&) = delete;
};
}  // namespace detail

/// Count `delta` (+1 when one starts, -1 when it ends) threads that work for
/// a solve beside its calling thread outside parallel_for, e.g. HiGHS's root
/// tasks (#163), in the live and peak counts with the spawned ones.
inline void count_busy_threads(int delta) {
    detail::add_live(delta);
}

/// Most threads parallel_for had spawned at once since the last reset, those
/// count_busy_threads reported included.
inline int peak_spawned_threads() {
    return detail::peak_spawned.load(std::memory_order_relaxed);
}

/// Restart the peak count from the threads running now.
inline void reset_peak_spawned_threads() {
    detail::peak_spawned.store(detail::live_spawned.load(std::memory_order_relaxed),
                               std::memory_order_relaxed);
}

/// Run f(i) for every i in [begin, end) on at most num_threads threads, the
/// calling thread among them: it spawns num_threads - 1, and each thread takes
/// the next index until none is left, so slow indices do not hold up others.
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

    std::atomic<int> next{begin};
    auto work = [&f, &next, end] {
        for (int i = next.fetch_add(1, std::memory_order_relaxed); i < end;
             i = next.fetch_add(1, std::memory_order_relaxed)) {
            f(i);
        }
    };
    std::vector<std::jthread> threads;
    threads.reserve(static_cast<size_t>(num_threads - 1));
    for (int t = 0; t < num_threads - 1; ++t) {
        threads.emplace_back([&work] {
            const detail::SpawnedThreadScope counted;
            work();
        });
    }
    work();
    // jthreads auto-join on destruction
}

/// Spawn tasks that auto-join on destruction (thin jthread wrapper).
class task_group {
    std::vector<std::jthread> threads_;

public:
    template <typename F>
    void run(F&& f) {
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
