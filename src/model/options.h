#pragma once

#include <algorithm>
#include <cctype>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace cptp {

/// Parse a boolean option value: true/false, 1/0 or on/off, in any case.
/// Throws std::invalid_argument naming the option on anything else, so a
/// typo or an unsupported spelling such as "yes" never reads as false.
inline bool parse_bool_option(std::string_view key, std::string_view value) {
    std::string lower(value);
    std::ranges::transform(lower, lower.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "true" || lower == "1" || lower == "on") {
        return true;
    }
    if (lower == "false" || lower == "0" || lower == "off") {
        return false;
    }
    throw std::invalid_argument(std::string(key) + " must be true/false, 1/0 or on/off, not '" +
                                std::string(value) + "'");
}

/// Parse a strictly positive int option value: the whole string must be a
/// decimal integer >= 1 that fits in an int. Throws std::invalid_argument
/// naming the option otherwise, so "0", "-5", "10x" or "" never reach code
/// that divides by the value.
inline int parse_positive_int_option(std::string_view key, std::string_view value) {
    int result = 0;
    const char* const first = value.data();
    const char* const end = first + value.size();
    auto [ptr, ec] = std::from_chars(first, end, result);
    if (ec != std::errc{} || ptr != end || result < 1) {
        throw std::invalid_argument(std::string(key) + " must be a positive integer, not '" +
                                    std::string(value) + "'");
    }
    return result;
}

}  // namespace cptp
