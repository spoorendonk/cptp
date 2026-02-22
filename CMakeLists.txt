cmake_minimum_required(VERSION 3.25)
project(rcspp-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/CompilerFlags.cmake)
include(cmake/FetchHiGHS.cmake)

# TBB — system package
find_package(TBB REQUIRED)

# ---------- separator library (solver-independent) ----------
add_library(rcspp_sep
    src/sep/sec_separator.cpp
    src/sep/rci_separator.cpp
    src/sep/multistar_separator.cpp
    src/sep/comb_separator.cpp
    src/sep/rglm_separator.cpp
)
target_include_directories(rcspp_sep PUBLIC src)
target_link_libraries(rcspp_sep PUBLIC TBB::tbb)
rcspp_set_compiler_flags(rcspp_sep)

# ---------- core library ----------
add_library(rcspp_core
    src/core/problem.cpp
    src/core/io.cpp
)
target_include_directories(rcspp_core PUBLIC src)
target_link_libraries(rcspp_core PUBLIC)
rcspp_set_compiler_flags(rcspp_core)

# ---------- model library (HiGHS integration) ----------
add_library(rcspp_model
    src/model/model.cpp
    src/model/highs_bridge.cpp
)
target_include_directories(rcspp_model PUBLIC src)
# HiGHS internal headers (mip/HighsUserSeparator.h) are PRIVATE in the highs target,
# so we need to add them explicitly for our cutpool callback integration.
target_include_directories(rcspp_model PRIVATE ${highs_SOURCE_DIR}/highs)
target_link_libraries(rcspp_model PUBLIC rcspp_core rcspp_sep highs)
rcspp_set_compiler_flags(rcspp_model)

# ---------- CLI ----------
add_executable(rcspp-solve src/cli/main.cpp)
target_link_libraries(rcspp-solve PRIVATE rcspp_model)

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

    add_executable(rcspp_tests
        tests/test_max_flow.cpp
        tests/test_separators.cpp
        tests/test_model.cpp
        tests/test_optimal.cpp
        tests/test_propagator.cpp
    )
    target_link_libraries(rcspp_tests PRIVATE rcspp_model Catch2::Catch2WithMain)
    catch_discover_tests(rcspp_tests)

    # Instance integration tests
    add_executable(test_instances tests/test_instances.cpp)
    target_link_libraries(test_instances PRIVATE rcspp_model Catch2::Catch2WithMain)
    rcspp_set_compiler_flags(test_instances)

    add_executable(test_formulation tests/test_formulation.cpp)
    target_link_libraries(test_formulation PRIVATE rcspp_model)
endif()
