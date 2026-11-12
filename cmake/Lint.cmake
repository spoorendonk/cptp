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
# `.venv/bin` is searched first, as the hooks do, so the tools the dev group
# pins there (`pyproject.toml`) win over the ones on PATH; NO_CACHE, so one put
# there after a first configure is picked up. The expected version is the
# group's one toolchain version, 22.1.8, read from `.clang-format-version`
# (which the pre-commit hook reads too). A missing tool, or one of another
# major version (it formats this tree differently, or reports other
# findings), is a configure warning that drops its tests, so a machine
# without them still builds; with the CI environment variable set, the wrong
# major is a configure error instead.
#
# Expects CPTP_LINT_UNITS (paths relative to the source root) to be set.

set(_cptp_lint_hint "${CMAKE_CURRENT_SOURCE_DIR}/.venv/bin")
find_program(CPTP_CLANG_FORMAT clang-format HINTS "${_cptp_lint_hint}" NO_CACHE)
find_program(CPTP_CLANG_TIDY clang-tidy HINTS "${_cptp_lint_hint}" NO_CACHE)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/.clang-format-version"
     CPTP_CLANG_VERSION LIMIT_COUNT 1)
string(STRIP "${CPTP_CLANG_VERSION}" CPTP_CLANG_VERSION)
string(REGEX MATCH "^([0-9]+)" _cptp_lint_major "${CPTP_CLANG_VERSION}")
if(NOT _cptp_lint_major)
    message(FATAL_ERROR ".clang-format-version does not begin with a version "
                        "(\"${CPTP_CLANG_VERSION}\")")
endif()

# Clears <var>, with a warning (an error when CI is set), when the tool's
# major version is not the pinned one.
function(_cptp_lint_require_pinned var)
    if(NOT ${var})
        return()
    endif()
    execute_process(COMMAND ${${var}} --version
        OUTPUT_VARIABLE _v ERROR_QUIET RESULT_VARIABLE _rc)
    string(REGEX MATCH "version ([0-9]+)\\." _m "${_v}")
    if(NOT _rc EQUAL 0 OR NOT CMAKE_MATCH_1 STREQUAL _cptp_lint_major)
        set(_level WARNING)
        if(DEFINED ENV{CI})
            set(_level FATAL_ERROR)
        endif()
        message(${_level} "${${var}} is not version ${CPTP_CLANG_VERSION} "
                          "(\"${_v}\"): its lint tests are not registered. "
                          "Run `uv sync` for the pinned copy in .venv/bin.")
        set(${var} "" PARENT_SCOPE)
    endif()
endfunction()
_cptp_lint_require_pinned(CPTP_CLANG_FORMAT)
_cptp_lint_require_pinned(CPTP_CLANG_TIDY)

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

# On macOS, Apple's compiler shims find the SDK (and with it libc++'s
# headers) through xcrun, so compile_commands.json carries no -isysroot when
# CMAKE_OSX_SYSROOT is empty. The pinned clang-tidy is a plain LLVM build that
# does not ask xcrun, and without a sysroot finds no standard headers at all.
# Hand it the SDK xcrun names. With CMAKE_OSX_SYSROOT set (a path or an SDK
# name), CMake already puts the resolved -isysroot in the database, so nothing
# is added. Other platforms are unaffected.
set(_cptp_tidy_extra_args "")
if(APPLE AND CPTP_CLANG_TIDY AND NOT CMAKE_OSX_SYSROOT)
    execute_process(COMMAND xcrun --show-sdk-path
        OUTPUT_VARIABLE _cptp_sdk OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET RESULT_VARIABLE _cptp_sdk_rc)
    if(_cptp_sdk_rc EQUAL 0 AND _cptp_sdk)
        list(APPEND _cptp_tidy_extra_args "--extra-arg=-isysroot${_cptp_sdk}")
    else()
        set(_level WARNING)
        if(DEFINED ENV{CI})
            set(_level FATAL_ERROR)
        endif()
        message(${_level} "No macOS SDK found (xcrun --show-sdk-path) and "
                          "CMAKE_OSX_SYSROOT is empty: clang-tidy would find "
                          "no standard headers, so its lint tests are not "
                          "registered.")
        set(CPTP_CLANG_TIDY "")
    endif()
endif()

if(CPTP_CLANG_TIDY)
    foreach(unit IN LISTS CPTP_LINT_UNITS)
        add_test(NAME "clang_tidy:${unit}"
            COMMAND ${CPTP_CLANG_TIDY} -p ${CMAKE_BINARY_DIR}
                    --quiet --warnings-as-errors=*
                    ${_cptp_tidy_extra_args}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${unit})
        set_tests_properties("clang_tidy:${unit}" PROPERTIES LABELS lint)
    endforeach()
else()
    message(WARNING "No usable clang-tidy: no clang_tidy: tests are "
                    "registered, so the push gate runs without lint.")
endif()
