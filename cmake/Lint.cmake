# Lint gates as ctest tests, so the pre-push local CI enforces them (#79).
#
#   clang_format        clang-format --dry-run -Werror over every tracked C++
#                       file, the same set the pre-commit hook checks. Needs a
#                       git checkout to list them, so a source export has none.
#   clang_tidy:<path>   clang-tidy on one translation unit, warnings as errors.
#                       One test per unit, named by its path from the
#                       repository root: pre-push then tidies only the units a
#                       push reaches (headers through a stub unit, via the
#                       ninja dependency log) and caches the passes.
#
# The units are every source the build compiles, so tidying them reaches every
# header; `.clang-tidy`'s HeaderFilterRegex reports on this repository's
# `src/`, `tests/` and `python/` and keeps out FetchContent's sources (HiGHS,
# Catch2, nanobind) under the build tree.
#
# `.venv/bin` is searched first, as the hooks do, so pinned tools there win
# over the ones on PATH; NO_CACHE, so one put there after a first configure is
# picked up. A missing tool, or one older than the 18 CLAUDE.md vouches for
# (older clang-format formats this tree differently), is a configure warning
# that drops its tests, not an error: a machine without them still builds.
#
# Expects CPTP_LINT_UNITS (paths relative to the source root) to be set.

set(_cptp_lint_hint "${CMAKE_CURRENT_SOURCE_DIR}/.venv/bin")
find_program(CPTP_CLANG_FORMAT clang-format HINTS "${_cptp_lint_hint}" NO_CACHE)
find_program(CPTP_CLANG_TIDY clang-tidy HINTS "${_cptp_lint_hint}" NO_CACHE)

# Clears <var> when the tool's major version is below 18, with a warning.
function(_cptp_lint_require_18 var)
    if(NOT ${var})
        return()
    endif()
    execute_process(COMMAND ${${var}} --version
        OUTPUT_VARIABLE _v ERROR_QUIET RESULT_VARIABLE _rc)
    string(REGEX MATCH "version ([0-9]+)\\." _m "${_v}")
    if(NOT _rc EQUAL 0 OR NOT CMAKE_MATCH_1 OR CMAKE_MATCH_1 LESS 18)
        message(WARNING "${${var}} is not version 18 or newer (\"${_v}\"): "
                        "its lint tests are not registered.")
        set(${var} "" PARENT_SCOPE)
    endif()
endfunction()
_cptp_lint_require_18(CPTP_CLANG_FORMAT)
_cptp_lint_require_18(CPTP_CLANG_TIDY)

if(CPTP_CLANG_FORMAT AND NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    message(STATUS "Not a git checkout: no clang_format test")
elseif(CPTP_CLANG_FORMAT)
    add_test(NAME clang_format
        COMMAND ${CMAKE_COMMAND}
            -DCLANG_FORMAT=${CPTP_CLANG_FORMAT}
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_clang_format.cmake)
    set_tests_properties(clang_format PROPERTIES LABELS lint)
else()
    message(WARNING "No usable clang-format: no clang_format test is "
                    "registered, so the push gate does not check formatting.")
endif()

# Units that include HiGHS's mip/HighsMipSolverData.h are tidied as C++20.
# Its HighsMipWorker.h resets a std::unique_ptr to a type it only forward
# declares in an inline destructor; under C++23, where unique_ptr is
# constexpr, clang rejects that as a hard error (gcc, which builds the
# project, accepts it), and no check option can silence an error.
set(CPTP_LINT_CXX20_UNITS src/model/model.cpp src/model/highs_bridge.cpp)

if(CPTP_CLANG_TIDY)
    foreach(unit IN LISTS CPTP_LINT_UNITS)
        set(_cptp_lint_std)
        if(unit IN_LIST CPTP_LINT_CXX20_UNITS)
            set(_cptp_lint_std --extra-arg=-std=c++20)
        endif()
        add_test(NAME "clang_tidy:${unit}"
            COMMAND ${CPTP_CLANG_TIDY} -p ${CMAKE_BINARY_DIR}
                    --quiet --warnings-as-errors=* ${_cptp_lint_std}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${unit})
        set_tests_properties("clang_tidy:${unit}" PROPERTIES LABELS lint)
    endforeach()
else()
    message(WARNING "No usable clang-tidy: no clang_tidy: tests are "
                    "registered, so the push gate runs without lint.")
endif()
