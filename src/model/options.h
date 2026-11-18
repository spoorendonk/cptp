#pragma once

#include <algorithm>
#include <cctype>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

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

/// Parse a root cut-loop limit (root_sepa_rounds, root_sepa_stall): "highs"
/// keeps HiGHS's own and is -1, "unlimited" lifts it and is 0, and a positive
/// integer N replaces it and is N (HighsUserHooks's encoding). Throws
/// std::invalid_argument naming the option on anything else.
inline int parse_root_limit_option(std::string_view key, std::string_view value) {
    if (value == "highs") {
        return -1;
    }
    if (value == "unlimited") {
        return 0;
    }
    try {
        return parse_positive_int_option(key, value);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument(std::string(key) +
                                    " must be highs, unlimited or a positive integer, not '" +
                                    std::string(value) + "'");
    }
}

/// The option value parse_root_limit_option reads as `limit`.
inline std::string root_limit_name(int limit) {
    if (limit < 0) {
        return "highs";
    }
    if (limit == 0) {
        return "unlimited";
    }
    return std::to_string(limit);
}

/// The limits of HiGHS's root cut loop on the main MIP (#150), each in
/// parse_root_limit_option's encoding: -1 keeps HiGHS's own (a cap of
/// 2*sqrt(log2 tree size) rounds, a stall limit of 3), 0 lifts it, N > 0
/// replaces it.
struct RootSeparationLimits {
    int rounds = -1;  // root_sepa_rounds: the round cap
    int stall = -1;   // root_sepa_stall: stalled rounds before the loop stops

    /// Parse `value` into the limit `key` names and return true, or return
    /// false when `key` is neither option.
    bool parse(std::string_view key, std::string_view value) {
        if (key == "root_sepa_rounds") {
            rounds = parse_root_limit_option(key, value);
            return true;
        }
        if (key == "root_sepa_stall") {
            stall = parse_root_limit_option(key, value);
            return true;
        }
        return false;
    }

    /// Append each limit that is not HiGHS's as "key=value".
    void append_nondefault(std::vector<std::string>& settings) const {
        if (rounds != -1) {
            settings.push_back("root_sepa_rounds=" + root_limit_name(rounds));
        }
        if (stall != -1) {
            settings.push_back("root_sepa_stall=" + root_limit_name(stall));
        }
    }
};

}  // namespace cptp
