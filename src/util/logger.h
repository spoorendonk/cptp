#pragma once

#include <format>
#include <iostream>
#include <mutex>
#include <ostream>
#include <string_view>

namespace rcspp {

class Logger {
 public:
    explicit Logger(std::ostream& out = std::cout) : out_(out) {}

    void log(std::string_view msg) {
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
    std::mutex mu_;
};

}  // namespace rcspp
