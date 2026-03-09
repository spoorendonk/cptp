#pragma once

#include <atomic>
#include <format>
#include <iostream>
#include <mutex>
#include <ostream>
#include <string_view>

namespace cptp {

class Logger {
 public:
  explicit Logger(std::ostream& out = std::cout) : out_(out) {}

  void set_enabled(bool enabled) {
    enabled_.store(enabled, std::memory_order_relaxed);
  }

  bool enabled() const { return enabled_.load(std::memory_order_relaxed); }

  void log(std::string_view msg) {
    if (!enabled()) return;
    std::lock_guard lock(mu_);
    out_ << msg;
    if (!msg.empty() && msg.back() != '\n') out_ << '\n';
    out_.flush();
  }

  template <typename... Args>
  void log(std::format_string<Args...> fmt, Args&&... args) {
    auto msg = std::format(fmt, std::forward<Args>(args)...);
    log(std::string_view{msg});
  }

 private:
  std::ostream& out_;
  std::atomic<bool> enabled_{true};
  std::mutex mu_;
};

}  // namespace cptp
