cmake_minimum_required(VERSION 3.25)
project(cptp-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/CompilerFlags.cmake)
include(cmake/FetchHiGHS.cmake)

# TBB — system package
find_package(TBB REQUIRED)

# ---------- separator library (solver-independent) ----------
add_library(cptp_sep
    src/sep/sec_separator.cpp
    src/sep/rci_separator.cpp
    src/sep/multistar_separator.cpp
    src/sep/comb_separator.cpp
)
target_include_directories(cptp_sep PUBLIC src)
target_link_libraries(cptp_sep PUBLIC TBB::tbb)
cptp_set_compiler_flags(cptp_sep)

# ---------- core library ----------
add_library(cptp_core
    src/core/problem.cpp
    src/core/io.cpp
)
target_include_directories(cptp_core PUBLIC src)
target_link_libraries(cptp_core PUBLIC)
cptp_set_compiler_flags(cptp_core)

# ---------- model library (HiGHS integration) ----------
add_library(cptp_model
    src/model/model.cpp
    src/model/highs_bridge.cpp
)
target_include_directories(cptp_model PUBLIC src)
# HiGHS internal headers (mip/HighsUserSeparator.h) are PRIVATE in the highs target,
# so we need to add them explicitly for our cutpool callback integration.
target_include_directories(cptp_model PRIVATE ${highs_SOURCE_DIR}/highs)
target_link_libraries(cptp_model PUBLIC cptp_core cptp_sep highs)
cptp_set_compiler_flags(cptp_model)

# ---------- CLI ----------
add_executable(cptp-solve src/cli/main.cpp)
target_link_libraries(cptp-solve PRIVATE cptp_model)

# ---------- Python bindings ----------
option(CPTP_BUILD_PYTHON "Build Python bindings" OFF)
if(CPTP_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# ---------- Tests ----------
option(CPTP_BUILD_TESTS "Build tests" ON)
if(CPTP_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.7.1
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)

    add_executable(cptp_tests
        tests/test_max_flow.cpp
        tests/test_separators.cpp
        tests/test_model.cpp
    )
    target_link_libraries(cptp_tests PRIVATE cptp_model Catch2::Catch2WithMain)
    catch_discover_tests(cptp_tests)

    # Benchmarks
    add_executable(bench_separators bench/bench_separators.cpp)
    target_link_libraries(bench_separators PRIVATE cptp_model Catch2::Catch2WithMain)
    cptp_set_compiler_flags(bench_separators)

    add_executable(test_formulation tests/test_formulation.cpp)
    target_link_libraries(test_formulation PRIVATE cptp_model)
endif()
