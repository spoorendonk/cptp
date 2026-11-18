#pragma once

#include <algorithm>
#include <string>
#include <thread>

namespace cptp::test {

// A test that solves a non-trivial instance under a time limit and requires
// optimality pins its threads option to at most kTestThreads and carries the
// [procs4] tag ([optimal] tests count too), for which CMake declares 4
// PROCESSORS. ctest -j then runs fewer of them side by side, and a
// time-limited solve stays clear of its limit (#160, #164). The default,
// (hw+1)/2, grows with the host; up to 8 logical CPUs it stays in force.
inline const std::string kTestThreads =
    std::to_string(std::min(4U, (std::max(1U, std::thread::hardware_concurrency()) + 1U) / 2U));

}  // namespace cptp::test
