#include "Highs.h"
#include "mip/HighsUserSeparator.h"

#include <algorithm>
#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace {

constexpr double kFeasTol = 1e-5;

// A random multi-dimensional knapsack with set-packing rows: hard enough for
// HiGHS to start sub-MIPs from inside sub-MIPs.
HighsLp nested_submip_model(int32_t num_col, int32_t num_row, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int32_t> weight(5, 60);
    std::uniform_int_distribution<int32_t> profit(10, 90);
    std::bernoulli_distribution in_row(0.35);

    HighsLp lp;
    lp.num_col_ = num_col;
    lp.num_row_ = num_row;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_.resize(num_col);
    lp.col_lower_.assign(num_col, 0.0);
    lp.col_upper_.assign(num_col, 1.0);
    lp.integrality_.assign(num_col, HighsVarType::kInteger);
    for (auto& cost : lp.col_cost_) {
        cost = profit(rng);
    }

    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.num_col_ = num_col;
    lp.a_matrix_.num_row_ = num_row;
    lp.a_matrix_.start_.assign(1, 0);
    lp.row_lower_.assign(num_row, -kHighsInf);
    lp.row_upper_.resize(num_row);
    for (int32_t row = 0; row < num_row; ++row) {
        double total = 0.0;
        for (int32_t col = 0; col < num_col; ++col) {
            if (in_row(rng)) {
                const double value = weight(rng);
                lp.a_matrix_.index_.push_back(col);
                lp.a_matrix_.value_.push_back(value);
                total += value;
            }
        }
        lp.row_upper_[row] = std::floor(total / 2.0);
        lp.a_matrix_.start_.push_back(static_cast<HighsInt>(lp.a_matrix_.index_.size()));
    }
    return lp;
}

bool satisfies_rows(const HighsLp& lp, const std::vector<double>& sol) {
    for (int32_t row = 0; row < lp.num_row_; ++row) {
        double activity = 0.0;
        for (HighsInt k = lp.a_matrix_.start_[row]; k < lp.a_matrix_.start_[row + 1]; ++k) {
            activity += lp.a_matrix_.value_[k] * sol[lp.a_matrix_.index_[k]];
        }
        if (activity > lp.row_upper_[row] + kFeasTol) {
            return false;
        }
    }
    return true;
}

// The knapsack above plus, for each of its columns x_j, a complement column
// z_j with x_j + z_j = 1: a doubleton equation, which a sub-MIP's presolve
// removes by substituting one column for the other.
HighsLp complement_model(int32_t num_col, int32_t num_row, uint32_t seed) {
    HighsLp lp = nested_submip_model(num_col, num_row, seed);
    std::mt19937 rng(seed + 1);
    std::uniform_int_distribution<int32_t> profit(1, 30);
    lp.num_col_ = 2 * num_col;
    lp.num_row_ = num_row + num_col;
    lp.a_matrix_.num_col_ = lp.num_col_;
    lp.a_matrix_.num_row_ = lp.num_row_;
    for (int32_t j = 0; j < num_col; ++j) {
        lp.col_cost_.push_back(profit(rng));
        lp.col_lower_.push_back(0.0);
        lp.col_upper_.push_back(1.0);
        lp.integrality_.push_back(HighsVarType::kInteger);
        lp.a_matrix_.index_.push_back(j);
        lp.a_matrix_.value_.push_back(1.0);
        lp.a_matrix_.index_.push_back(num_col + j);
        lp.a_matrix_.value_.push_back(1.0);
        lp.a_matrix_.start_.push_back(static_cast<HighsInt>(lp.a_matrix_.index_.size()));
        lp.row_lower_.push_back(1.0);
        lp.row_upper_.push_back(1.0);
    }
    return lp;
}

}  // namespace

TEST_CASE("Sub-MIP cut translation: a cut column without an image drops the cut",
          "[highs_patch][submip]") {
    using Image = HighsUserSeparator::ColumnImage;
    const std::vector<Image> images{
        {0, 1.0, 0.0},
        {HighsUserSeparator::kImageNone, 0.0, 0.0},
        {HighsUserSeparator::kImageConstant, 0.0, 1.0},
    };
    std::vector<HighsInt> idx{0, 1};
    std::vector<double> val{1.0, 1.0};
    double rhs = 1.0;
    CHECK_FALSE(HighsUserSeparator::toSolveSpace(images, idx, val, rhs));

    // An index outside the images is refused too.
    std::vector<HighsInt> far_idx{3};
    std::vector<double> far_val{1.0};
    double far_rhs = 1.0;
    CHECK_FALSE(HighsUserSeparator::toSolveSpace(images, far_idx, far_val, far_rhs));
}

TEST_CASE("Sub-MIP cut translation: fixed columns fold, shared columns merge",
          "[highs_patch][submip]") {
    using Image = HighsUserSeparator::ColumnImage;
    // x0 = y0; x1 = 1 - y0 (a doubleton equation's partner); x2 fixed at 1;
    // x3 = 2 * y1 + 0.5; x4 = -0.5 * y2 + 1.
    const std::vector<Image> images{
        {0, 1.0, 0.0}, {0, -1.0, 1.0}, {HighsUserSeparator::kImageConstant, 0.0, 1.0},
        {1, 2.0, 0.5}, {2, -0.5, 1.0},
    };

    SECTION("coefficients on one column sum, and a cancelled one is dropped") {
        // x0 + x1 + 3 x2 + x3 <= 6  ->  2 y1 <= 6 - 1 - 3 - 0.5
        std::vector<HighsInt> idx{0, 1, 2, 3};
        std::vector<double> val{1.0, 1.0, 3.0, 1.0};
        double rhs = 6.0;
        REQUIRE(HighsUserSeparator::toSolveSpace(images, idx, val, rhs));
        CHECK(idx == std::vector<HighsInt>{1});
        CHECK(val == std::vector<double>{2.0});
        CHECK(rhs == 1.5);
    }
    SECTION("distinct columns map in order, with scale and constant") {
        // 2 x4 + x0 - x1 <= 1  ->  -y2 + 2 y0 <= 1 - 2 + 1
        std::vector<HighsInt> idx{4, 0, 1};
        std::vector<double> val{2.0, 1.0, -1.0};
        double rhs = 1.0;
        REQUIRE(HighsUserSeparator::toSolveSpace(images, idx, val, rhs));
        CHECK(idx == std::vector<HighsInt>{0, 2});
        CHECK(val == std::vector<double>{2.0, -1.0});
        CHECK(rhs == 0.0);
    }
}

TEST_CASE("Sub-MIP cut translation: an image undoPrimal disagrees with is cleared",
          "[highs_patch][submip]") {
    using Image = HighsUserSeparator::ColumnImage;
    std::vector<Image> images{
        {0, 1.0, 0.0},
        {0, -1.0, 1.0},
        {HighsUserSeparator::kImageConstant, 0.0, 1.0},
        {1, 1.0, 0.0},
        {HighsUserSeparator::kImageNone, 0.0, 0.0},
    };
    const std::vector<double> sol{0.25, 0.5};
    // x1 should be 0.75 and x2 1.0; the expansion says otherwise for both.
    const std::vector<double> expanded{0.25, 0.5, 0.0, 0.5, 0.3};
    CHECK(HighsUserSeparator::clearMismatchedImages(images, sol, expanded, 1e-9) == 2);
    CHECK(images[0].col == 0);
    CHECK(images[1].col == HighsUserSeparator::kImageNone);
    CHECK(images[2].col == HighsUserSeparator::kImageNone);
    CHECK(images[3].col == 1);
    CHECK(images[4].col == HighsUserSeparator::kImageNone);
}

TEST_CASE("Sub-MIP cut translation: images of fixed and substituted columns match undoPrimal",
          "[highs_patch][submip]") {
    const HighsLp lp = complement_model(60, 20, 11);
    const auto num_top = static_cast<size_t>(lp.num_col_);

    std::atomic<int64_t> rounds{0};
    std::atomic<int64_t> substituted{0};
    std::atomic<int64_t> fixed{0};
    std::atomic<int64_t> bad_images{0};
    std::atomic<int64_t> bad_translations{0};

    HighsUserSeparator::setCallback([&](const HighsLpRelaxation& relaxation, HighsCutPool&,
                                        const HighsMipSolver& mipsolver) {
        if (!mipsolver.submip) {
            return;
        }
        rounds++;
        const auto& reduced = relaxation.getSolution().col_value;
        const auto expanded = HighsUserSeparator::toOriginalSpace(mipsolver, reduced);
        auto images = HighsUserSeparator::columnImages(mipsolver);
        if (expanded.size() != num_top || images.size() != num_top) {
            bad_images++;
            return;
        }
        // Every image gives the column's expanded value.
        bad_images += HighsUserSeparator::clearMismatchedImages(images, reduced, expanded, 1e-9);

        // A cut over every column with an image, coefficient j + 1.
        std::vector<HighsInt> idx;
        std::vector<double> val;
        double activity = 0.0;
        for (size_t j = 0; j < num_top; ++j) {
            const auto& image = images[j];
            if (image.col == HighsUserSeparator::kImageNone) {
                continue;
            }
            if (image.col == HighsUserSeparator::kImageConstant) {
                fixed++;
            } else if (image.scale != 1.0 || image.constant != 0.0) {
                substituted++;  // a doubleton equation's substituted column
            }
            idx.push_back(static_cast<HighsInt>(j));
            val.push_back(static_cast<double>(j + 1));
            activity += static_cast<double>(j + 1) * expanded[j];
        }

        // Translated, that cut has the same slack at the reduced solution.
        const double rhs = activity + 0.5;
        double red_rhs = rhs;
        if (!HighsUserSeparator::toSolveSpace(images, idx, val, red_rhs)) {
            bad_translations++;
            return;
        }
        double red_activity = 0.0;
        for (size_t k = 0; k < idx.size(); ++k) {
            red_activity += val[k] * reduced[static_cast<size_t>(idx[k])];
        }
        if (std::abs((red_rhs - red_activity) - (rhs - activity)) > 1e-6 * (1.0 + std::abs(rhs))) {
            bad_translations++;
        }
    });

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("mip_heuristic_effort", 0.8);
    highs.setOptionValue("time_limit", 60.0);
    const HighsStatus pass_status = highs.passModel(lp);
    const HighsStatus run_status =
        pass_status == HighsStatus::kOk ? highs.run() : HighsStatus::kError;
    // Clear before any REQUIRE: the hooks capture this frame's locals.
    HighsUserSeparator::clearCallback();
    REQUIRE(pass_status == HighsStatus::kOk);
    REQUIRE(run_status != HighsStatus::kError);

    // The test proves nothing unless sub-MIPs ran, with presolve both fixing
    // columns and substituting them through the complement rows.
    REQUIRE(rounds.load() > 0);
    REQUIRE(substituted.load() > 0);
    REQUIRE(fixed.load() > 0);
    CHECK(bad_images.load() == 0);
    CHECK(bad_translations.load() == 0);
}

TEST_CASE("Nested sub-MIP: callbacks see the top-level column space", "[highs_patch][submip]") {
    const HighsLp lp = nested_submip_model(60, 20, 7);

    std::atomic<HighsInt> deepest_level{0};
    std::atomic<int64_t> nested_rounds{0};
    std::atomic<int64_t> bad_expansions{0};
    std::atomic<int64_t> bad_columns{0};
    std::atomic<int64_t> checks{0};
    std::atomic<int64_t> bad_checks{0};

    HighsUserSeparator::setCallback([&](const HighsLpRelaxation& relaxation, HighsCutPool&,
                                        const HighsMipSolver& mipsolver) {
        HighsInt seen = deepest_level.load();
        while (seen < mipsolver.submip_level &&
               !deepest_level.compare_exchange_weak(seen, mipsolver.submip_level)) {
        }
        if (mipsolver.submip_level < 2) {
            return;
        }
        nested_rounds++;

        const auto& reduced = relaxation.getSolution().col_value;
        const auto expanded = HighsUserSeparator::toOriginalSpace(mipsolver, reduced);
        if (std::cmp_not_equal(expanded.size(), lp.num_col_) || !satisfies_rows(lp, expanded)) {
            bad_expansions++;
            return;
        }
        const auto cols =
            HighsUserSeparator::originalColumns(mipsolver, static_cast<HighsInt>(reduced.size()));
        for (size_t r = 0; r < reduced.size(); ++r) {
            if (cols[r] < 0 || cols[r] >= lp.num_col_ ||
                std::abs(expanded[static_cast<size_t>(cols[r])] - reduced[r]) > kFeasTol) {
                bad_columns++;
            }
        }
    });
    HighsUserSeparator::setFeasibilityCheck([&](const std::vector<double>& sol) {
        checks++;
        if (std::cmp_not_equal(sol.size(), lp.num_col_) || !satisfies_rows(lp, sol)) {
            bad_checks++;
        }
        return true;
    });

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("mip_heuristic_effort", 0.8);
    highs.setOptionValue("time_limit", 60.0);
    const HighsStatus pass_status = highs.passModel(lp);
    const HighsStatus run_status =
        pass_status == HighsStatus::kOk ? highs.run() : HighsStatus::kError;
    // Clear before any REQUIRE: the hooks capture this frame's locals, and
    // a later test on this thread must not run them.
    HighsUserSeparator::clearCallback();
    REQUIRE(pass_status == HighsStatus::kOk);
    REQUIRE(run_status != HighsStatus::kError);

    // Without a level-2 sub-MIP the test would prove nothing.
    REQUIRE(deepest_level.load() >= 2);
    REQUIRE(nested_rounds.load() > 0);
    CHECK(bad_expansions.load() == 0);
    CHECK(bad_columns.load() == 0);
    REQUIRE(checks.load() > 0);
    CHECK(bad_checks.load() == 0);
}
