#include "Highs.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsUserSeparator.h"

#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <functional>
#include <thread>
#include <vector>

// HiGHS's root tasks (symmetry detection, the analytic centre) reported to the
// solve's root task observer while they run beside the root (#163).

namespace {

struct RootTasksRun {
    int reports = 0;          // tasks reported running on another thread
    int max_running = 0;      // most reported running at once
    int pending_at_sep = -1;  // reported pending at the first root cut-loop call
    int running_at_sep = -1;  // reported running then
    int pending_after = -1;   // reported pending after the solve
    int running_after = -1;   // reported running after the solve
    bool separated = false;   // the separator ran in the root cut loop
    HighsStatus status = HighsStatus::kError;
};

// Poll until done() holds, for at most `limit`; whether it did.
bool wait_until(const std::function<bool()>& done, std::chrono::milliseconds limit) {
    const auto deadline = std::chrono::steady_clock::now() + limit;
    while (!done()) {
        if (std::chrono::steady_clock::now() >= deadline) {
            return false;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return true;
}

// max x + y over x, y in [0, 2] with an integer z that nothing constrains.
// The lazy check rejects every solution until the separator has run in the
// main MIP's root cut loop, so the root LP does not close the solve before it.
// With hold, a task reported running waits (bounded) in the observer until
// that call, and the call waits (bounded) for such a report: the call then runs
// beside a task HiGHS has not joined yet, as root separation does beside a
// slow one.
// HiGHS records rootlpsol right after the loop, so an empty one marks a call in
// it.
RootTasksRun run_root(int threads, bool hold) {
    HighsLp lp;
    lp.num_col_ = 3;
    lp.num_row_ = 0;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = {1.0, 1.0, 0.0};
    lp.col_lower_ = {0.0, 0.0, 0.0};
    lp.col_upper_ = {2.0, 2.0, 1.0};
    lp.integrality_ = {HighsVarType::kContinuous, HighsVarType::kContinuous,
                       HighsVarType::kInteger};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 3;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0, 0, 0};

    RootTasksRun run;
    std::atomic<int> reports{0};
    std::atomic<int> pending{0};
    std::atomic<int> running{0};
    std::atomic<int> max_running{0};
    std::atomic<bool> separated{false};
    HighsUserSeparator::setRootTaskObserver([&, hold](int pending_delta, int running_delta) {
        pending += pending_delta;
        const int now = running.fetch_add(running_delta) + running_delta;
        if (running_delta > 0) {
            ++reports;
            int seen = max_running.load();
            while (now > seen && !max_running.compare_exchange_weak(seen, now)) {
            }
            if (hold) {
                wait_until([&] { return separated.load(); }, std::chrono::seconds(10));
            }
        }
    });
    HighsUserSeparator::setCallback(
        [&, hold](const HighsLpRelaxation&, HighsCutPool&, const HighsMipSolver& mipsolver) {
            if (mipsolver.submip || !mipsolver.mipdata_->rootlpsol.empty() || separated.load()) {
                return;
            }
            if (hold) {
                wait_until([&] { return running.load() > 0; }, std::chrono::seconds(2));
            }
            run.pending_at_sep = pending.load();
            run.running_at_sep = running.load();
            separated = true;
        });
    HighsUserSeparator::setFeasibilityCheck(
        [&separated](const std::vector<double>&) { return separated.load(); });

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("threads", threads);
    highs.setOptionValue("time_limit", 60.0);
    // An earlier solve on this thread may have left a scheduler of another
    // size, and HiGHS refuses another `threads` against it.
    Highs::resetThreadScheduler(true);
    if (highs.passModel(lp) == HighsStatus::kOk) {
        run.status = highs.run();
    }
    // Before any REQUIRE: the hooks capture this frame's state, and a later
    // test on this thread must not run them.
    HighsUserSeparator::clearCallback();
    run.reports = reports.load();
    run.max_running = max_running.load();
    run.pending_after = pending.load();
    run.running_after = running.load();
    run.separated = separated.load();
    return run;
}

}  // namespace

TEST_CASE("Root tasks: the observer is pending until clearCallback resets it",
          "[highs_patch][root_tasks]") {
    HighsUserSeparator::clearCallback();
    REQUIRE(HighsUserHooks::capture() == nullptr);

    // An observer alone makes a solve capture the hooks.
    HighsUserSeparator::setRootTaskObserver([](int, int) {});
    const auto hooks = HighsUserHooks::capture();
    // Before any REQUIRE, so a failure leaves no observer pending.
    HighsUserSeparator::clearCallback();
    REQUIRE(hooks != nullptr);
    CHECK(static_cast<bool>(hooks->root_task_observer));

    CHECK(!HighsUserHooks::pending().root_task_observer);
    REQUIRE(HighsUserHooks::capture() == nullptr);
}

// threads = 2 gives HiGHS one worker beside the solving thread, which almost
// always takes a spawned root task at once: root separation then sees it
// running. HiGHS's scheduler does not promise it (in about 2 runs in 100 the
// worker misses the task, and the solving thread runs it itself at the join,
// rightly not reported running), so the test tries a few solves for one where
// it does. [procs4] reserves processors under ctest -j, so a loaded machine
// does not keep HiGHS's worker from the task.
TEST_CASE("Root tasks: a task on another thread is reported while root separation runs",
          "[highs_patch][root_tasks][procs4]") {
    RootTasksRun run;
    for (int attempt = 0; attempt < 5 && run.running_at_sep < 1; ++attempt) {
        run = run_root(2, /*hold=*/true);
    }
    REQUIRE(run.status == HighsStatus::kOk);
    REQUIRE(run.separated);
    CHECK(run.reports >= 1);
    CHECK(run.pending_at_sep >= run.running_at_sep);
    CHECK(run.running_at_sep >= 1);
    CHECK(run.max_running <= 1);  // HiGHS's one worker
    CHECK(run.pending_after == 0);
    CHECK(run.running_after == 0);
}

// threads = 1: HiGHS has no worker, so the solving thread runs each root task
// itself when it joins it, after the root cut loop, and no other thread is kept
// busy. The task is still pending at the cut loop: spawned, not yet started, as
// one a worker has yet to wake for is, so the callbacks leave it its thread.
TEST_CASE("Root tasks: a task not yet started is pending, one the solving thread runs not running",
          "[highs_patch][root_tasks]") {
    const auto run = run_root(1, /*hold=*/false);
    REQUIRE(run.status == HighsStatus::kOk);
    REQUIRE(run.separated);
    CHECK(run.reports == 0);
    CHECK(run.pending_at_sep >= 1);
    CHECK(run.running_at_sep == 0);
    CHECK(run.pending_after == 0);
    CHECK(run.running_after == 0);
}
