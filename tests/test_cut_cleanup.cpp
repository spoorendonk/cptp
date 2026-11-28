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
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x);
    REQUIRE(lost);
    CHECK(lost == 0.0);
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{1.0, 1.0});
    CHECK(rhs == 1.0);
}

TEST_CASE("Cut cleanup: a 1e-17 that should be zero drops", "[cut_cleanup]") {
    // A GLM's 1 - 2*q/Q for q = Q/2 can come out as ~1e-17 instead of 0.
    std::vector<int32_t> idx{0, 1, 2};
    std::vector<double> val{1.0, 1e-17, 1.0};
    double rhs = 1.0;
    const std::vector<double> x{1.0, 0.5, 0.5};
    REQUIRE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x));
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{1.0, 1.0});
    CHECK(rhs == 1.0);  // at lower bound 0
}

TEST_CASE("Cut cleanup: a small positive coefficient moves to the rhs at its lower bound",
          "[cut_cleanup]") {
    // Column 1 in [1, 5]: the lower bound relaxes a positive term, the
    // upper one would cut off feasible points.
    const std::vector<double> lower{0.0, 1.0};
    const std::vector<double> upper{1.0, 5.0};
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, 1e-10};
    double rhs = 1.0;
    const std::vector<double> x{1.0, 3.0};
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, lower, upper, x);
    REQUIRE(lost);
    CHECK(idx == std::vector<int32_t>{0});
    CHECK(rhs == 1.0 - 1e-10 * 1.0);
    CHECK(lost == 1e-10 * (3.0 - 1.0));
}

TEST_CASE("Cut cleanup: a small negative coefficient on a binary relaxes the rhs",
          "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1, 2};
    std::vector<double> val{2.0, -1e-10, 2.0};
    double rhs = 2.0;
    const std::vector<double> x{1.0, 0.25, 0.5};
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x);
    REQUIRE(lost);
    CHECK(idx == std::vector<int32_t>{0, 2});
    CHECK(val == std::vector<double>{2.0, 2.0});
    // -a * ub = 1e-10 * 1: every point the cut allowed it still allows.
    CHECK(rhs == 2.0 - (-1e-10 * 1.0));
    CHECK(lost == -1e-10 * (0.25 - 1.0));
}

TEST_CASE("Cut cleanup: a coefficient above the threshold stays", "[cut_cleanup]") {
    // An RCI coefficient 2*q/Q with q/Q = 5e-7 is small but meant.
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, 1e-6};
    double rhs = 0.5;
    const std::vector<double> x{1.0, 1.0};
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x);
    REQUIRE(lost);
    CHECK(lost == 0.0);
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
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, upper, x));
}

TEST_CASE("Cut cleanup: a cut left empty is skipped", "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{0.0, 0.0};
    double rhs = -1.0;
    const std::vector<double> x{1.0, 1.0, 1.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x));
}

TEST_CASE("Cut cleanup: a cut whose coefficients are all noise is skipped", "[cut_cleanup]") {
    // Relative to max|a| nothing would drop; the threshold is absolute.
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1e-16, -1e-16};
    double rhs = 0.0;
    const std::vector<double> x{1.0, 0.0, 0.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x));
}

TEST_CASE("Cut cleanup: a +-1 term stays beside a huge coefficient", "[cut_cleanup]") {
    // RCI's 2*q_i/Q_r is ~1e10 for Q_r = 2e-10 * q_i; a threshold relative to
    // max|a| would drop the -1 edge terms and lose violation of order 1.
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{-1.0, 1e10};
    double rhs = 0.0;
    const std::vector<double> x{0.5, 0.5};
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x);
    CHECK(lost == 0.0);
    CHECK(idx == std::vector<int32_t>{0, 1});
    CHECK(rhs == 0.0);
}

TEST_CASE("Cut cleanup: a NaN coefficient skips the cut", "[cut_cleanup]") {
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, std::numeric_limits<double>::quiet_NaN()};
    double rhs = 1.0;
    const std::vector<double> x{1.0, 1.0, 0.0};
    CHECK_FALSE(drop_near_zero_coefficients(idx, val, rhs, kLower, kUpper, x));
}

TEST_CASE("Cut cleanup: the violation lost at the LP point is returned", "[cut_cleanup]") {
    // x0 - 1e-10 x1 <= 1 with x1 <= 1e9: dropping x1 relaxes the rhs by 0.1,
    // and at x1 = 0 the point's violation falls by as much.
    const std::vector<double> lower{0.0, 0.0};
    const std::vector<double> upper{1e9, 1e9};
    std::vector<int32_t> idx{0, 1};
    std::vector<double> val{1.0, -1e-10};
    double rhs = 1.0;
    const std::vector<double> x{1.15, 0.0};
    const auto lost = drop_near_zero_coefficients(idx, val, rhs, lower, upper, x);
    REQUIRE(lost);
    CHECK(idx == std::vector<int32_t>{0});
    CHECK(rhs == 1.0 - (-1e-10 * 1e9));
    CHECK(lost == -1e-10 * (0.0 - 1e9));
}
