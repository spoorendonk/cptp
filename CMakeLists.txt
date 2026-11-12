cmake_minimum_required(VERSION 3.25)
project(cptp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/CompilerFlags.cmake)

# ============================================================
# cptp — unified solver library
# Core + algorithms + HiGHS integration
# ============================================================
include(cmake/FetchHiGHS.cmake)

add_library(cptp
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
    src/heuristic/primal_heuristic.cpp
    src/model/model.cpp
    src/model/highs_bridge.cpp
    src/model/hyperplane.cpp
)
target_include_directories(cptp PUBLIC src)
# HiGHS internal headers (mip/HighsUserSeparator.h) are PRIVATE in the highs target,
# so we need to add them explicitly for our cutpool callback integration.
target_include_directories(cptp PRIVATE ${highs_SOURCE_DIR}/highs)
target_link_libraries(cptp PUBLIC highs)
cptp_set_compiler_flags(cptp)
# Embed the project version (single source of truth: project(... VERSION ...)).
# PRIVATE so it reaches the library's own sources (model.cpp preamble) without
# colliding with the generated version.h the CLI includes.
target_compile_definitions(cptp PRIVATE CPTP_VERSION="${PROJECT_VERSION}")

# ---------- CLI ----------
# Regenerate version.h on every build so the embedded git hash tracks HEAD
# without a reconfigure. configure_file is content-idempotent, so this avoids
# spurious relinks when the hash hasn't changed.
set(CPTP_VERSION_HEADER ${CMAKE_CURRENT_BINARY_DIR}/generated/version.h)
add_custom_target(cptp_version ALL
    BYPRODUCTS ${CPTP_VERSION_HEADER}
    COMMAND ${CMAKE_COMMAND}
        -DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DCPTP_VERSION=${PROJECT_VERSION}
        -DVERSION_IN=${CMAKE_CURRENT_SOURCE_DIR}/src/cli/version.h.in
        -DVERSION_OUT=${CPTP_VERSION_HEADER}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_version.cmake
    COMMENT "Generating cptp version.h"
    VERBATIM)

add_executable(cptp-solve src/cli/main.cpp)
add_dependencies(cptp-solve cptp_version)
target_link_libraries(cptp-solve PRIVATE cptp)
target_include_directories(cptp-solve PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)

# ---------- Tests ----------
option(CPTP_BUILD_TESTS "Build tests" ON)

# ---------- Python bindings ----------
option(CPTP_BUILD_PYTHON "Build Python bindings" OFF)

# The binding tests run in ctest (python_tests) when an interpreter with pytest
# and numpy is found, `.venv/bin` first, and the bindings are then built here
# too. Skipped with a warning otherwise: the C++ build must not depend on a
# Python setup. NO_CACHE, so a .venv made after a first configure is picked up.
# Top level only, and never inside a wheel build (SKBUILD).
set(CPTP_PYTHON_TESTS OFF)
if(CPTP_BUILD_TESTS AND PROJECT_IS_TOP_LEVEL AND NOT SKBUILD)
    find_program(CPTP_TEST_PYTHON NAMES python3 python
        HINTS "${CMAKE_CURRENT_SOURCE_DIR}/.venv/bin" NAMES_PER_DIR NO_CACHE)
    if(CPTP_TEST_PYTHON)
        # 3.10+: the package uses `X | None` annotations.
        execute_process(COMMAND ${CPTP_TEST_PYTHON} -c
            "import sys, numpy, pytest; sys.exit(sys.version_info < (3, 10))"
            RESULT_VARIABLE CPTP_TEST_PYTHON_OK OUTPUT_QUIET ERROR_QUIET)
    endif()
    if(CPTP_TEST_PYTHON AND CPTP_TEST_PYTHON_OK EQUAL 0)
        # The interpreter python/ builds the bindings for.
        set(Python_EXECUTABLE ${CPTP_TEST_PYTHON})
        find_package(Python COMPONENTS Interpreter Development.Module)
    endif()
    if(Python_FOUND)
        set(CPTP_PYTHON_TESTS ON)
    else()
        message(WARNING "No Python 3.10+ with pytest, numpy and its development "
                        "headers found: python_tests is not registered, so the "
                        "gate skips the binding tests. Install with: uv sync")
    endif()
endif()

if(CPTP_BUILD_PYTHON OR CPTP_PYTHON_TESTS)
    add_subdirectory(python)
endif()
if(CPTP_BUILD_TESTS)
    # FetchHiGHS sets BUILD_TESTING=OFF in cache for third-party builds.
    # Re-enable testing in this directory scope without forcing global cache.
    set(BUILD_TESTING ON)
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
    add_executable(cptp_tests
        tests/test_model.cpp
        tests/test_multicap.cpp
        tests/test_optimal.cpp
        tests/test_path_cut_validity.cpp
        tests/test_propagator.cpp
        tests/test_spi.cpp
        tests/test_thread_local.cpp
        tests/test_ubs_cutoff.cpp
    )
    target_include_directories(cptp_tests PRIVATE ${highs_SOURCE_DIR}/highs)
    target_link_libraries(cptp_tests PRIVATE cptp Catch2::Catch2WithMain)
    cptp_set_compiler_flags(cptp_tests)
    catch_discover_tests(cptp_tests
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    # Extended unit tests (kept separate to control runtime impact).
    add_executable(cptp_tests_extra
        tests/test_max_flow.cpp
        tests/test_oracle.cpp
        tests/test_separators.cpp
    )
    target_link_libraries(cptp_tests_extra PRIVATE cptp Catch2::Catch2WithMain)
    cptp_set_compiler_flags(cptp_tests_extra)
    catch_discover_tests(cptp_tests_extra
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    # Instance integration tests
    add_executable(test_instances tests/test_instances.cpp)
    target_link_libraries(test_instances PRIVATE cptp Catch2::Catch2WithMain)
    cptp_set_compiler_flags(test_instances)

    add_executable(test_formulation tests/test_formulation.cpp)
    target_link_libraries(test_formulation PRIVATE cptp)
    cptp_set_compiler_flags(test_formulation)

    # Heuristic-stage UB regression gate over benchmark corpora.
    # Intentionally not auto-discovered in ctest due long runtime.
    add_executable(test_heuristic_stage_ub tests/test_heuristic_stage_ub.cpp)
    target_link_libraries(test_heuristic_stage_ub PRIVATE cptp Catch2::Catch2WithMain)
    cptp_set_compiler_flags(test_heuristic_stage_ub)

    if(CPTP_PYTHON_TESTS)
        add_test(NAME python_tests
            COMMAND ${CMAKE_COMMAND}
                -DPYTHON=${CPTP_TEST_PYTHON}
                -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DSTAGE_DIR=${CMAKE_CURRENT_BINARY_DIR}/python_stage
                -DEXTENSION=$<TARGET_FILE:cptp_python>
                -DSOLVE=$<TARGET_FILE:cptp-solve>
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_pytest.cmake)
        # The committed type stubs match what this build's extension exports.
        add_test(NAME python_stub
            COMMAND ${CMAKE_COMMAND}
                -DGENERATED=${CPTP_STUB}
                -DCOMMITTED=${CPTP_STUB_SOURCE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_stub.cmake)
    endif()
endif()

# Lint gates. Top level only: they read our compile_commands.json, which only a
# top-level build writes, and a consumer's ctest is not ours to lint.
if(CPTP_BUILD_TESTS AND PROJECT_IS_TOP_LEVEL)
    # Every translation unit the build compiles; tidying them reaches every
    # header. See cmake/Lint.cmake.
    set(CPTP_LINT_UNITS)
    foreach(target IN ITEMS cptp cptp-solve cptp_tests cptp_tests_extra
            test_instances test_formulation test_heuristic_stage_ub cptp_python)
        if(NOT TARGET ${target})
            continue()
        endif()
        get_target_property(target_sources ${target} SOURCES)
        get_target_property(target_dir ${target} SOURCE_DIR)
        foreach(src IN LISTS target_sources)
            if(src MATCHES "\\.cpp$")
                cmake_path(ABSOLUTE_PATH src BASE_DIRECTORY ${target_dir})
                cmake_path(RELATIVE_PATH src BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
                list(APPEND CPTP_LINT_UNITS ${src})
            endif()
        endforeach()
    endforeach()
    include(cmake/Lint.cmake)
endif()
