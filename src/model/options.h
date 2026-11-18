#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
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

/// The on/off switches of the cut families, each at its default.
struct CutFamilyFlags {
    bool enable_sec = true;
    bool enable_rci = true;
    bool enable_multistar = true;
    bool enable_comb = false;
    bool enable_rglm = false;
    bool enable_spi = false;
    // RCI and RGLM also separate on the surrogate row of all bounded rows
    // (#112), next to the capacity rows. Off until measured (#113).
    bool capacity_cut_surrogate = false;
    // Multistar also separates on each target's largest depot-target min-cut
    // side, not only the Gomory-Hu tree's smallest (#168). Off until measured.
    bool multistar_largest_side = false;

    /// Parse `value` into the switch `key` names and return true, or return
    /// false when `key` names none. Throws as parse_bool_option does.
    bool parse(std::string_view key, std::string_view value) {
        const auto all = flags();
        const auto* it =
            std::ranges::find_if(all, [key](const auto& entry) { return entry.first == key; });
        if (it == all.end()) {
            return false;
        }
        this->*(it->second) = parse_bool_option(key, value);
        return true;
    }

    /// Append each switch that is not at its default as "key=on" or "key=off".
    void append_nondefault(std::vector<std::string>& settings) const {
        const CutFamilyFlags defaults;
        for (const auto& [name, flag] : flags()) {
            if (this->*flag != defaults.*flag) {
                settings.push_back(std::string(name) + (this->*flag ? "=on" : "=off"));
            }
        }
    }

private:
    using Flag = bool CutFamilyFlags::*;
    static std::array<std::pair<std::string_view, Flag>, 8> flags() {
        return {{
            {"enable_sec", &CutFamilyFlags::enable_sec},
            {"enable_rci", &CutFamilyFlags::enable_rci},
            {"enable_multistar", &CutFamilyFlags::enable_multistar},
            {"enable_comb", &CutFamilyFlags::enable_comb},
            {"enable_rglm", &CutFamilyFlags::enable_rglm},
            {"enable_spi", &CutFamilyFlags::enable_spi},
            {"capacity_cut_surrogate", &CutFamilyFlags::capacity_cut_surrogate},
            {"multistar_largest_side", &CutFamilyFlags::multistar_largest_side},
        }};
    }
};

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
