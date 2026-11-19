#include "core/io.h"
#include "model/model.h"
#include "parallel/parallel.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <latch>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>

using Catch::Matchers::WithinAbs;

namespace {

// The distinct threads parallel_for ran the body on, and whether every index
// ran exactly once.
struct Observed {
    std::set<std::thread::id> threads;
    bool each_once = true;
};

Observed run_parallel_for(int range, int num_threads) {
    std::mutex mutex;
    std::vector<int> hits(static_cast<size_t>(range), 0);
    Observed observed;
    cptp::parallel::parallel_for(
        0, range,
        [&](int i) {
            const std::scoped_lock lock(mutex);
            observed.threads.insert(std::this_thread::get_id());
            ++hits[static_cast<size_t>(i)];
        },
        num_threads);
    for (const int h : hits) {
        observed.each_once = observed.each_once && h == 1;
    }
    return observed;
}

}  // namespace

// The calling thread is one of the num_threads, but it may find every index
// already taken by the spawned ones, so the test does not require it ran one.
TEST_CASE("parallel_for: runs each index once on at most num_threads threads", "[parallel]") {
    const int num_threads = GENERATE(2, 3, 5);
    INFO("num_threads = " << num_threads);
    cptp::parallel::reset_peak_spawned_threads();
    const auto observed = run_parallel_for(64, num_threads);
    REQUIRE(observed.each_once);
    REQUIRE(observed.threads.size() <= static_cast<size_t>(num_threads));
    REQUIRE(cptp::parallel::peak_spawned_threads() <= num_threads - 1);
}

// No thread count falls back to hardware_concurrency (#159).
TEST_CASE("parallel_for: num_threads <= 1 runs on the calling thread", "[parallel]") {
    const int num_threads = GENERATE(-1, 0, 1);
    INFO("num_threads = " << num_threads);
    cptp::parallel::reset_peak_spawned_threads();
    const auto observed = run_parallel_for(64, num_threads);
    REQUIRE(observed.each_once);
    REQUIRE(observed.threads == std::set{std::this_thread::get_id()});
    REQUIRE(cptp::parallel::peak_spawned_threads() == 0);
}

// The latch holds every thread inside the body until all four are there, so
// the three spawned ones run at once and the peak must see exactly them.
TEST_CASE("parallel_for: peak counts spawned threads running at once", "[parallel]") {
    cptp::parallel::reset_peak_spawned_threads();
    std::latch all_in(4);
    cptp::parallel::parallel_for(0, 4, [&](int) { all_in.arrive_and_wait(); }, 4);
    REQUIRE(cptp::parallel::peak_spawned_threads() == 3);
    cptp::parallel::reset_peak_spawned_threads();
    REQUIRE(cptp::parallel::peak_spawned_threads() == 0);
}

// A solve keeps its parallel_for threads within `threads` (#159): the
// separators, RCI's candidate loop, reduced-cost fixing (fix-to-one too), the
// LP-guided heuristic and the all-pairs labeling share the budget, the calling
// thread among them, with the HiGHS root tasks (symmetry detection, the
// analytic centre) running beside the root's callbacks (#163), so at most
// threads - 1 threads run beside the calling one at once. parallel=off keeps
// HiGHS's parallel lock inactive, so every callback may use the whole budget
// those tasks leave. Not counted: HiGHS's tree-search workers and stage-1's
// forward/backward labeling thread.
TEST_CASE("Model: a solve's parallel_for spawns at most threads - 1 threads",
          "[parallel][model][slow][procs4]") {
    const int threads = GENERATE(1, 3);
    const std::string parallel = GENERATE("off", "on");
    INFO("threads = " << threads << ", parallel = " << parallel);
    cptp::Model model;
    model.set_problem(cptp::io::load("benchmarks/instances/spprclib/B-n45-k6-54.sppcc"));
    const cptp::SolverOptions opts = {
        {"output_flag", "false"}, {"time_limit", "60"},      {"threads", std::to_string(threads)},
        {"parallel", parallel},   {"rc_fixing", "adaptive"}, {"rc_fixing_to_one", "true"}};
    cptp::parallel::reset_peak_spawned_threads();
    const auto r = model.solve(opts);
    REQUIRE(r.is_optimal());
    REQUIRE_THAT(r.objective, WithinAbs(-74278.0, 1.0));
    REQUIRE(cptp::parallel::peak_spawned_threads() <= threads - 1);
}
