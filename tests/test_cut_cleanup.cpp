#include "sep/cut_cleanup.h"

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <vector>

using cptp::sep::drop_near_zero_coefficients;

namespace {

// Three binary columns, unless a test says otherwise.
const std::vector<double> kLower{0.0, 0.0, 0.0};
const std::vector<double> kUpper{1.0, 1.0, 1.0};

}  // namespace

TEST_CASE("Cut cleanup: an exact zero drops and leaves the rhs", "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1, 2};
    std::vector<double> val{1.0, 0.0, 1.0};
    double rhs = 1.0;
    const std::vector<double> x{1.0, 0.5, 0.5};
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x, 0.1));
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{1.0, 1.0});
    CHECK(rhs == 1.0);
}

TEST_CASE("Cut cleanup: a 1e-17 that should be zero drops at its lower bound", "[cut_cleanup]") {
    // A GLM's 1 - 2*q/Q for q = Q/2 can come out as ~1e-17 instead of 0.
    // Column 1 has lower bound 1 here, so dropping it moves 1e-17 to the rhs.
    const std::vector<double> lower{0.0, 1.0, 0.0};
    std::vector<int32_t> idx{0, 1, 2};
    std::vector<double> val{1.0, 1e-17, 1.0};
    double rhs = 1.0;
    const std::vector<double> x{1.0, 1.0, 0.5};
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, lower, kUpper, x, 0.1));
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{1.0, 1.0});
    CHECK(rhs == 1.0 - 1e-17);
}

TEST_CASE("Cut cleanup: a small negative coefficient on a binary relaxes the rhs",
          "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1, 2};
    std::vector<double> val{2.0, -1e-10, 2.0};
    double rhs = 2.0;
    const std::vector<double> x{1.0, 1.0, 0.5};
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x, 0.1));
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{2.0, 2.0});
    // -a * ub = 1e-10 * 1: every point the cut allowed it still allows.
    CHECK(rhs == 2.0 + 1e-10);
}

TEST_CASE("Cut cleanup: a coefficient above the threshold stays", "[cut_cleanup]") {
    // An RCI coefficient 2*q/Q with q/Q = 5e-7 is small but meant.
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, 1e-6};
    double rhs = 0.5;
    const std::vector<double> x{1.0, 1.0};
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x, 0.1));
    CHECK(idx == std::vector<int32_t>{0, 1});
    CHECK(val == std::vector<double>{1.0, 1e-6});
    CHECK(rhs == 0.5);
}

TEST_CASE("Cut cleanup: a term to drop with an infinite bound skips the cut", "[cut_cleanup]") {
    const std::vector<double> upper{1.0, std::numeric_limits<double>::infinity(), 1.0};
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, -1e-12};
    double rhs = 0.5;
    const std::vector<double> x{1.0, 0.0, 0.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, upper, x, 0.1));
}

TEST_CASE("Cut cleanup: a cut left empty is skipped", "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{0.0, 0.0};
    double rhs = -1.0;
    const std::vector<double> x{1.0, 1.0, 1.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x, 0.1));
}

TEST_CASE("Cut cleanup: a cut that no longer cuts off the LP point is skipped", "[cut_cleanup]") {
    // x0 - 1e-10 x1 <= 1 with x1 <= 1e9: dropping x1 relaxes the rhs by 0.1,
    // and the LP point's violation falls from 0.15 to 0.05.
    const std::vector<double> upper{1e9, 1e9};
    const std::vector<double> lower{0.0, 0.0};
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, -1e-10};
    double rhs = 1.0;
    const std::vector<double> x{1.15, 0.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, lower, upper, x, 0.1));

    // The same cut with a smaller minimum violation stays.
    idx = {0, 1};
    val = {1.0, -1e-10};
    rhs = 1.0;
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, lower, upper, x, 0.01));
    CHECK(idx == std::vector<int32_t>{0});
    CHECK(rhs == 1.0 - (-1e-10 * 1e9));
}
