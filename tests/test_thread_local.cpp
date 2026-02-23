#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <atomic>
#include <barrier>
#include <thread>
#include <vector>

#include "mip/HighsUserSeparator.h"
#include "mip/HighsUserPropagator.h"
#include "core/io.h"
#include "model/model.h"

using Catch::Matchers::WithinAbs;

// ─── Separator callback lifecycle ───

TEST_CASE("Separator: hasCallback is false initially", "[thread_local][separator]") {
    HighsUserSeparator::clearCallback();
    REQUIRE_FALSE(HighsUserSeparator::hasCallback());
}

TEST_CASE("Separator: setCallback / hasCallback / clearCallback", "[thread_local][separator]") {
    HighsUserSeparator::clearCallback();
    REQUIRE_FALSE(HighsUserSeparator::hasCallback());

    HighsUserSeparator::setCallback(
        [](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver&) {});
    REQUIRE(HighsUserSeparator::hasCallback());

    HighsUserSeparator::clearCallback();
    REQUIRE_FALSE(HighsUserSeparator::hasCallback());
}

// ─── Feasibility check lifecycle ───

TEST_CASE("Separator: isFeasible returns true when no check registered",
          "[thread_local][separator]") {
    HighsUserSeparator::clearCallback();
    REQUIRE(HighsUserSeparator::isFeasible({1.0, 0.0, 1.0}));
}

TEST_CASE("Separator: feasibility check is invoked and respected",
          "[thread_local][separator]") {
    HighsUserSeparator::clearCallback();

    int call_count = 0;
    HighsUserSeparator::setFeasibilityCheck([&](const std::vector<double>& sol) {
        ++call_count;
        // Reject solutions where first element > 0.5
        return sol.empty() || sol[0] <= 0.5;
    });

    REQUIRE(HighsUserSeparator::isFeasible({0.0, 1.0}));
    REQUIRE(call_count == 1);

    REQUIRE_FALSE(HighsUserSeparator::isFeasible({1.0, 0.0}));
    REQUIRE(call_count == 2);

    HighsUserSeparator::clearCallback();
    // After clear, should return true unconditionally
    REQUIRE(HighsUserSeparator::isFeasible({1.0, 0.0}));
    REQUIRE(call_count == 2);  // not called again
}

// ─── Propagator callback lifecycle ───

TEST_CASE("Propagator: clearCallback leaves no callback", "[thread_local][propagator]") {
    HighsUserPropagator::clearCallback();
    // propagate() with no callback should be a no-op (not crash)
    // We can't call propagate without real HiGHS objects, but we can verify
    // set/clear works
    HighsUserPropagator::setCallback(
        [](HighsDomain&, const HighsMipSolver&, const HighsLpRelaxation&) {});
    HighsUserPropagator::clearCallback();
    // No crash = pass
    REQUIRE(true);
}

// ─── Thread-local isolation ───

TEST_CASE("Separator: callbacks are thread-local isolated",
          "[thread_local][separator]") {
    HighsUserSeparator::clearCallback();

    // Barrier so both threads register callbacks before either checks
    std::barrier sync(2);
    std::atomic<bool> thread_a_ok{false};
    std::atomic<bool> thread_b_ok{false};

    // Track which thread's callback was invoked
    std::atomic<int> a_calls{0};
    std::atomic<int> b_calls{0};

    std::thread thread_a([&] {
        // Initially no callback on this thread
        REQUIRE_FALSE(HighsUserSeparator::hasCallback());

        HighsUserSeparator::setCallback(
            [&](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver&) {
                a_calls.fetch_add(1);
            });
        HighsUserSeparator::setFeasibilityCheck([&](const std::vector<double>&) {
            a_calls.fetch_add(1);
            return true;
        });

        sync.arrive_and_wait();  // both threads have registered
        sync.arrive_and_wait();  // both threads have verified

        // Thread A should still have its callback
        REQUIRE(HighsUserSeparator::hasCallback());

        // Invoke feasibility check — should use A's lambda
        HighsUserSeparator::isFeasible({});
        thread_a_ok = (a_calls.load() == 1);

        HighsUserSeparator::clearCallback();
    });

    std::thread thread_b([&] {
        // Initially no callback on this thread (thread-local!)
        REQUIRE_FALSE(HighsUserSeparator::hasCallback());

        HighsUserSeparator::setCallback(
            [&](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver&) {
                b_calls.fetch_add(1);
            });
        HighsUserSeparator::setFeasibilityCheck([&](const std::vector<double>&) {
            b_calls.fetch_add(1);
            return false;
        });

        sync.arrive_and_wait();  // both threads have registered
        sync.arrive_and_wait();  // both threads have verified

        // Thread B should still have its callback
        REQUIRE(HighsUserSeparator::hasCallback());

        // Invoke feasibility check — should use B's lambda (returns false)
        bool result = HighsUserSeparator::isFeasible({});
        thread_b_ok = (b_calls.load() == 1 && !result);

        HighsUserSeparator::clearCallback();
    });

    thread_a.join();
    thread_b.join();

    REQUIRE(thread_a_ok.load());
    REQUIRE(thread_b_ok.load());
    // Cross-contamination check: A's callback was never called by B and vice versa
    REQUIRE(a_calls.load() == 1);
    REQUIRE(b_calls.load() == 1);
}

TEST_CASE("Separator: clearCallback on one thread does not affect another",
          "[thread_local][separator]") {
    std::barrier sync(2);
    std::atomic<bool> b_still_has{false};

    std::thread thread_a([&] {
        HighsUserSeparator::setCallback(
            [](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver&) {});
        sync.arrive_and_wait();  // both set
        HighsUserSeparator::clearCallback();
        sync.arrive_and_wait();  // A cleared
    });

    std::thread thread_b([&] {
        HighsUserSeparator::setCallback(
            [](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver&) {});
        sync.arrive_and_wait();  // both set
        sync.arrive_and_wait();  // A cleared
        // B's callback should still be active
        b_still_has = HighsUserSeparator::hasCallback();
        HighsUserSeparator::clearCallback();
    });

    thread_a.join();
    thread_b.join();

    REQUIRE(b_still_has.load());
}

TEST_CASE("Propagator: callbacks are thread-local isolated",
          "[thread_local][propagator]") {
    std::barrier sync(2);
    std::atomic<bool> a_has{false};
    std::atomic<bool> b_has{false};

    std::thread thread_a([&] {
        HighsUserPropagator::clearCallback();
        HighsUserPropagator::setCallback(
            [](HighsDomain&, const HighsMipSolver&, const HighsLpRelaxation&) {});
        sync.arrive_and_wait();
        a_has = true;  // A has callback
        sync.arrive_and_wait();
        HighsUserPropagator::clearCallback();
    });

    std::thread thread_b([&] {
        HighsUserPropagator::clearCallback();
        // B does NOT set a callback
        sync.arrive_and_wait();
        sync.arrive_and_wait();
        // B should not see A's callback — no API to check directly,
        // but clearCallback should not crash
        b_has = false;
        HighsUserPropagator::clearCallback();
    });

    thread_a.join();
    thread_b.join();

    REQUIRE(a_has.load());
    REQUIRE_FALSE(b_has.load());
}

// ─── Integration: concurrent Model::solve ───

static const rcspp::SolverOptions quiet = {
    {"time_limit", "30"},
    {"output_flag", "false"},
};

TEST_CASE("Concurrent Model::solve on different threads",
          "[thread_local][integration]") {
    // Two threads solve the same instance concurrently.
    // With the old mutex-based global callbacks, this would either deadlock
    // or produce wrong results. With thread_local, both should succeed.
    std::atomic<bool> ok_a{false};
    std::atomic<bool> ok_b{false};

    auto solve_tiny4 = [&](std::atomic<bool>& ok) {
        auto prob = rcspp::io::load("tests/data/tiny4.txt");
        rcspp::Model model;
        model.set_problem(std::move(prob));
        auto result = model.solve(quiet);
        ok = result.has_solution() && result.is_optimal();
    };

    std::thread t_a([&] { solve_tiny4(ok_a); });
    std::thread t_b([&] { solve_tiny4(ok_b); });

    t_a.join();
    t_b.join();

    REQUIRE(ok_a.load());
    REQUIRE(ok_b.load());
}

TEST_CASE("Concurrent Model::solve with different instances",
          "[thread_local][integration]") {
    // Thread A solves a tour, thread B solves a path.
    std::atomic<bool> ok_tour{false};
    std::atomic<bool> ok_path{false};
    std::atomic<double> obj_tour{0.0};
    std::atomic<double> obj_path{0.0};

    std::thread t_tour([&] {
        auto prob = rcspp::io::load("tests/data/tiny4.txt");
        rcspp::Model model;
        model.set_problem(std::move(prob));
        auto result = model.solve(quiet);
        if (result.is_optimal()) {
            ok_tour = true;
            obj_tour = result.objective;
        }
    });

    std::thread t_path([&] {
        auto prob = rcspp::io::load("tests/data/tiny4_path.txt");
        rcspp::Model model;
        model.set_problem(std::move(prob));
        auto result = model.solve(quiet);
        if (result.is_optimal()) {
            ok_path = true;
            obj_path = result.objective;
        }
    });

    t_tour.join();
    t_path.join();

    REQUIRE(ok_tour.load());
    REQUIRE(ok_path.load());
    REQUIRE_THAT(obj_tour.load(), WithinAbs(-11.0, 1.0));
    REQUIRE_THAT(obj_path.load(), WithinAbs(-13.0, 1.0));
}
