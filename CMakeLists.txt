cmake_minimum_required(VERSION 3.25)
project(rcspp-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/CompilerFlags.cmake)

# TBB — system package
find_package(TBB REQUIRED)

# ============================================================
# rcspp — unified solver library
# Core + algorithms + HiGHS integration
# ============================================================
include(cmake/FetchHiGHS.cmake)

add_library(rcspp
    src/core/problem.cpp
    src/core/io.cpp
    src/sep/sec_separator.cpp
    src/sep/rci_separator.cpp
    src/sep/multistar_separator.cpp
    src/sep/comb_separator.cpp
    src/sep/rglm_separator.cpp
    src/sep/spi_separator.cpp
    src/sep/separation_oracle.cpp
    src/preprocess/bound_propagator.cpp
    src/model/model.cpp
    src/model/highs_bridge.cpp
)
target_include_directories(rcspp PUBLIC src)
# HiGHS internal headers (mip/HighsUserSeparator.h) are PRIVATE in the highs target,
# so we need to add them explicitly for our cutpool callback integration.
target_include_directories(rcspp PRIVATE ${highs_SOURCE_DIR}/highs)
target_link_libraries(rcspp PUBLIC TBB::tbb highs)
rcspp_set_compiler_flags(rcspp)

# ---------- CLI ----------
add_executable(rcspp-solve src/cli/main.cpp)
target_link_libraries(rcspp-solve PRIVATE rcspp)

# ---------- Python bindings ----------
option(RCSPP_BUILD_PYTHON "Build Python bindings" OFF)
if(RCSPP_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# ---------- Tests ----------
option(RCSPP_BUILD_TESTS "Build tests" ON)
if(RCSPP_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.7.1
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)

    # HiGHS integration tests (only solver-dependent tests)
    add_executable(rcspp_tests
        tests/test_model.cpp
        tests/test_ng_labeling.cpp
        tests/test_optimal.cpp
        tests/test_path_cut_validity.cpp
        tests/test_propagator.cpp
        tests/test_spi.cpp
        tests/test_thread_local.cpp
        tests/test_ubs_cutoff.cpp
    )
    target_include_directories(rcspp_tests PRIVATE ${highs_SOURCE_DIR}/highs)
    target_link_libraries(rcspp_tests PRIVATE rcspp Catch2::Catch2WithMain)
    catch_discover_tests(rcspp_tests)

    # Instance integration tests
    add_executable(test_instances tests/test_instances.cpp)
    target_link_libraries(test_instances PRIVATE rcspp Catch2::Catch2WithMain)
    rcspp_set_compiler_flags(test_instances)

    add_executable(test_formulation tests/test_formulation.cpp)
    target_link_libraries(test_formulation PRIVATE rcspp)
endif()
