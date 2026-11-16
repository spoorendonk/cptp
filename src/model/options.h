#pragma once

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

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

}  // namespace cptp
