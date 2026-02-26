#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>

#include "core/io.h"
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
    opts.max_labels_per_node = 50;

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
