#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>

#include "core/io.h"
#include "model/async_incumbent.h"
#include "preprocess/edge_elimination.h"
#include "preprocess/ng_labeling.h"
#include "preprocess/shared_bounds.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("ng DSSR with ng_size=1 matches 2-cycle labeling baseline",
          "[labeling][ng]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");

    auto fwd_old = rcspp::preprocess::forward_labeling(prob, prob.source());
    auto bwd_old = rcspp::preprocess::backward_labeling(prob, prob.target());

    rcspp::preprocess::ng::DssrOptions opts;
    opts.initial_ng_size = 1;
    opts.max_ng_size = 1;
    opts.dssr_iterations = 1;

    auto bounds = rcspp::preprocess::ng::compute_bounds(
        prob, prob.source(), prob.target(), opts);

    REQUIRE(bounds.fwd.size() == fwd_old.size());
    REQUIRE(bounds.bwd.size() == bwd_old.size());
    for (size_t i = 0; i < fwd_old.size(); ++i) {
        REQUIRE_THAT(bounds.fwd[i], WithinAbs(fwd_old[i], 1e-9));
        REQUIRE_THAT(bounds.bwd[i], WithinAbs(bwd_old[i], 1e-9));
    }
}

TEST_CASE("shared bounds store tightens monotonically", "[labeling][shared_bounds]") {
    rcspp::preprocess::BoundSnapshot init{
        .fwd = {1.0, 2.0, 3.0},
        .bwd = {1.5, 2.5, 3.5},
        .correction = 0.0,
        .ng_size = 4,
        .version = 0,
    };

    rcspp::preprocess::SharedBoundsStore store(init);

    SECTION("weaker update does not change snapshot") {
        bool changed = store.publish_tightening(
            {0.5, 1.0, 2.0},
            {1.0, 2.0, 3.0},
            3);
        REQUIRE_FALSE(changed);
        auto snap = store.snapshot();
        REQUIRE(snap.version == 0);
        REQUIRE(snap.ng_size == 4);
        REQUIRE_THAT(snap.fwd[0], WithinAbs(1.0, 1e-12));
        REQUIRE_THAT(snap.bwd[2], WithinAbs(3.5, 1e-12));
    }

    SECTION("tighter update increases version and bounds") {
        bool changed = store.publish_tightening(
            {1.2, 2.0, 3.4},
            {1.5, 2.9, 3.7},
            8);
        REQUIRE(changed);
        auto snap = store.snapshot();
        REQUIRE(snap.version == 1);
        REQUIRE(snap.ng_size == 8);
        REQUIRE_THAT(snap.fwd[0], WithinAbs(1.2, 1e-12));
        REQUIRE_THAT(snap.fwd[1], WithinAbs(2.0, 1e-12));
        REQUIRE_THAT(snap.fwd[2], WithinAbs(3.4, 1e-12));
        REQUIRE_THAT(snap.bwd[1], WithinAbs(2.9, 1e-12));
    }
}

TEST_CASE("ng DSSR reports finite elementary path bound on tiny path",
          "[labeling][ng]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");
    rcspp::preprocess::ng::DssrOptions opts;
    opts.initial_ng_size = 4;
    opts.max_ng_size = 8;
    opts.dssr_iterations = 4;

    auto bounds = rcspp::preprocess::ng::compute_bounds(
        prob, prob.source(), prob.target(), opts);

    REQUIRE(bounds.fwd.size() == static_cast<size_t>(prob.num_nodes()));
    REQUIRE(bounds.bwd.size() == static_cast<size_t>(prob.num_nodes()));
    REQUIRE(bounds.fwd[prob.target()] < std::numeric_limits<double>::infinity());
    REQUIRE(bounds.elementary_path_found);
    REQUIRE(bounds.elementary_path_cost < std::numeric_limits<double>::infinity());
    REQUIRE_FALSE(bounds.elementary_path.empty());
    REQUIRE(bounds.elementary_path.front() == prob.source());
    REQUIRE(bounds.elementary_path.back() == prob.target());
}

TEST_CASE("ng DSSR accepts ng_initial_size=0",
          "[labeling][ng]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");

    rcspp::preprocess::ng::DssrOptions opts0;
    opts0.initial_ng_size = 0;
    opts0.max_ng_size = 1;
    opts0.dssr_iterations = 1;

    rcspp::preprocess::ng::DssrOptions opts1;
    opts1.initial_ng_size = 1;
    opts1.max_ng_size = 1;
    opts1.dssr_iterations = 1;

    auto b0 = rcspp::preprocess::ng::compute_bounds(
        prob, prob.source(), prob.target(), opts0);
    auto b1 = rcspp::preprocess::ng::compute_bounds(
        prob, prob.source(), prob.target(), opts1);

    REQUIRE(b0.fwd.size() == b1.fwd.size());
    REQUIRE(b0.bwd.size() == b1.bwd.size());
    for (size_t i = 0; i < b0.fwd.size(); ++i) {
        REQUIRE_THAT(b0.fwd[i], WithinAbs(b1.fwd[i], 1e-9));
        REQUIRE_THAT(b0.bwd[i], WithinAbs(b1.bwd[i], 1e-9));
    }
}

TEST_CASE("ng DSSR with zero-init still yields finite target bound",
          "[labeling][ng]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");
    rcspp::preprocess::ng::DssrOptions opts;
    opts.initial_ng_size = 0;
    opts.max_ng_size = 3;
    opts.dssr_iterations = 4;

    auto bounds = rcspp::preprocess::ng::compute_bounds(
        prob, prob.source(), prob.target(), opts);

    REQUIRE(bounds.fwd.size() == static_cast<size_t>(prob.num_nodes()));
    REQUIRE(bounds.bwd.size() == static_cast<size_t>(prob.num_nodes()));
    REQUIRE(bounds.fwd[prob.target()] < std::numeric_limits<double>::infinity());
    REQUIRE(bounds.ng_size >= 1);
    REQUIRE(bounds.ng_size <= 3);
}

TEST_CASE("ng label update marks predecessor bit (prevents 3-cycle return to root)",
          "[labeling][ng]") {
    auto prob = rcspp::io::load("tests/data/tiny4_path.txt");
    const int32_t n = prob.num_nodes();
    REQUIRE(n >= 3);

    std::vector<std::vector<int32_t>> ng_neighbors(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) {
        ng_neighbors[static_cast<size_t>(i)] = {0, 1, 2};
    }

    auto run = rcspp::preprocess::ng::detail::run_labeling(
        prob, 0, prob.edge_costs(), prob.profits(), ng_neighbors, true);

    // With predecessor-marking update, node 0 cannot be revisited via 0->1->2->0.
    REQUIRE_THAT(run.best_cost[0], WithinAbs(-prob.profit(0), 1e-9));
}

TEST_CASE("async incumbent store keeps best objective", "[labeling][async_incumbent]") {
    rcspp::model::AsyncIncumbentStore store;
    REQUIRE(store.publish_if_better({1.0, 0.0, 1.0}, -10.0));
    REQUIRE_FALSE(store.publish_if_better({0.0, 1.0, 0.0}, -5.0));
    REQUIRE(store.publish_if_better({0.0, 1.0, 1.0}, -12.0));
    auto snap = store.snapshot();
    REQUIRE_THAT(snap.objective, WithinAbs(-12.0, 1e-12));
    REQUIRE(snap.version == 2);
    REQUIRE(snap.col_values.size() == 3);
}
