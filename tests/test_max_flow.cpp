#include "core/digraph.h"
#include "core/dinitz.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <thread>

using Catch::Matchers::WithinAbs;

TEST_CASE("Dinitz max-flow on simple graph", "[maxflow]") {
    // Simple graph: 0 -> 1 -> 3, 0 -> 2 -> 3
    // Capacities: 0->1: 10, 1->3: 5, 0->2: 8, 2->3: 7
    // Max flow from 0 to 3 should be 12 (5 via 1, 7 via 2)
    cptp::digraph_builder builder(4);
    builder.add_arc(0, 1, 10.0);
    builder.add_arc(1, 0, 0.0);
    builder.add_arc(1, 3, 5.0);
    builder.add_arc(3, 1, 0.0);
    builder.add_arc(0, 2, 8.0);
    builder.add_arc(2, 0, 0.0);
    builder.add_arc(2, 3, 7.0);
    builder.add_arc(3, 2, 0.0);

    auto [graph, capacity] = builder.build();

    cptp::dinitz alg(graph, capacity, 0, 3);
    alg.run();

    REQUIRE_THAT(alg.flow_value(), WithinAbs(12.0, 1e-9));
}

TEST_CASE("Dinitz max-flow on classic example", "[maxflow]") {
    // Classic Ford-Fulkerson example:
    //   0 -> 1 (16), 0 -> 2 (13)
    //   1 -> 2 (4),  1 -> 3 (12)
    //   2 -> 1 (10), 2 -> 4 (14)
    //   3 -> 2 (9),  3 -> 5 (20)
    //   4 -> 3 (7),  4 -> 5 (4)
    // Max flow from 0 to 5 = 23
    cptp::digraph_builder builder(6);
    builder.add_arc(0, 1, 16.0);
    builder.add_arc(1, 0, 0.0);
    builder.add_arc(0, 2, 13.0);
    builder.add_arc(2, 0, 0.0);
    builder.add_arc(1, 2, 4.0);
    builder.add_arc(2, 1, 10.0);
    builder.add_arc(1, 3, 12.0);
    builder.add_arc(3, 1, 0.0);
    builder.add_arc(2, 4, 14.0);
    builder.add_arc(4, 2, 0.0);
    builder.add_arc(3, 2, 9.0);
    builder.add_arc(2, 3, 0.0);
    builder.add_arc(3, 5, 20.0);
    builder.add_arc(5, 3, 0.0);
    builder.add_arc(4, 3, 7.0);
    builder.add_arc(3, 4, 0.0);
    builder.add_arc(4, 5, 4.0);
    builder.add_arc(5, 4, 0.0);

    auto [graph, capacity] = builder.build();

    cptp::dinitz alg(graph, capacity, 0, 5);
    alg.run();

    REQUIRE_THAT(alg.flow_value(), WithinAbs(23.0, 1e-9));
}

TEST_CASE("Dinitz min-cut source side", "[maxflow]") {
    // 0 -> 1 (1), 0 -> 2 (1), 1 -> 3 (1), 2 -> 3 (1)
    // Max flow = 2
    cptp::digraph_builder builder(4);
    builder.add_arc(0, 1, 1.0);
    builder.add_arc(1, 0, 0.0);
    builder.add_arc(0, 2, 1.0);
    builder.add_arc(2, 0, 0.0);
    builder.add_arc(1, 3, 1.0);
    builder.add_arc(3, 1, 0.0);
    builder.add_arc(2, 3, 1.0);
    builder.add_arc(3, 2, 0.0);

    auto [graph, capacity] = builder.build();

    cptp::dinitz alg(graph, capacity, 0, 3);
    alg.run();

    REQUIRE_THAT(alg.flow_value(), WithinAbs(2.0, 1e-9));

    // Source side should include node 0, target side should include node 3
    REQUIRE(alg.on_source_side(0));
    REQUIRE(!alg.on_source_side(3));
}

TEST_CASE("Dinitz on symmetric graph (like SEC support)", "[maxflow]") {
    // Simulate undirected graph: 0-1 (cap 1), 0-2 (cap 1), 1-3 (cap 0.5), 2-3
    // (cap 0.5) Each undirected edge becomes two directed arcs
    cptp::digraph_builder builder(4);
    builder.add_arc(0, 1, 1.0);
    builder.add_arc(1, 0, 1.0);
    builder.add_arc(0, 2, 1.0);
    builder.add_arc(2, 0, 1.0);
    builder.add_arc(1, 3, 0.5);
    builder.add_arc(3, 1, 0.5);
    builder.add_arc(2, 3, 0.5);
    builder.add_arc(3, 2, 0.5);

    auto [graph, capacity] = builder.build();

    cptp::dinitz alg(graph, capacity, 0, 3);
    alg.run();

    // Max flow should be 2 * 0.5 = 1.0 (bottleneck at edges to node 3)
    REQUIRE_THAT(alg.flow_value(), WithinAbs(1.0, 1e-9));
}

TEST_CASE("Dinitz on a long path runs on a default-stack thread", "[maxflow]") {
    // A path-like support puts t n-1 BFS levels from s. The recursive search
    // this replaced took a stack frame per level and overflowed a default
    // std::thread stack (512 KB on macOS, where HiGHS's workers run it) at a
    // few thousand nodes (#161); 100000 levels is ~22 MB of those frames.
    constexpr int32_t kNodes = 100000;
    constexpr int32_t kNarrow = kNodes / 2;
    cptp::digraph_builder builder(kNodes);
    for (int32_t i = 0; i + 1 < kNodes; ++i) {
        const double cap = i == kNarrow ? 0.25 : 1.0;
        builder.add_arc(i, i + 1, cap);
        builder.add_arc(i + 1, i, cap);
    }
    auto [graph, capacity] = builder.build();

    double flow = -1.0;
    bool narrow_on_source_side = false;
    bool after_narrow_on_source_side = true;
    std::thread worker([&] {
        cptp::dinitz alg(graph, capacity, 0, kNodes - 1);
        alg.run();
        flow = alg.flow_value();
        narrow_on_source_side = alg.on_source_side(kNarrow);
        after_narrow_on_source_side = alg.on_source_side(kNarrow + 1);
    });
    worker.join();

    REQUIRE_THAT(flow, WithinAbs(0.25, 1e-9));
    CHECK(narrow_on_source_side);
    CHECK_FALSE(after_narrow_on_source_side);
}
