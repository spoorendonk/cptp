#pragma once

#include <algorithm>
#include <thread>
#include <vector>

namespace cptp::parallel {

/// Split [begin, end) across jthreads. num_threads=0 means hardware_concurrency.
template <typename F>
void parallel_for(int begin, int end, F&& f, int num_threads = 0) {
  if (begin >= end) return;
  if (num_threads <= 0)
    num_threads =
        static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
  const int range = end - begin;
  num_threads = std::min(num_threads, range);

  std::vector<std::jthread> threads;
  threads.reserve(static_cast<size_t>(num_threads));
  const int chunk = range / num_threads;
  int remainder = range % num_threads;

  int lo = begin;
  for (int t = 0; t < num_threads; ++t) {
    int hi = lo + chunk + (t < remainder ? 1 : 0);
    threads.emplace_back([&f, lo, hi] {
      for (int i = lo; i < hi; ++i) f(i);
    });
    lo = hi;
  }
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
      if (t.joinable()) t.join();
    }
    threads_.clear();
  }

  ~task_group() {
    for (auto& t : threads_) {
      if (t.joinable()) t.join();
    }
  }
};

}  // namespace cptp::parallel
