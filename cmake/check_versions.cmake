# The version_sync ctest test (CMakeLists.txt): CMake's project(VERSION) is
# the one source of cptp's version, and pyproject.toml, CITATION.cff and the
# README's software citation and version-DOI sentence each say the same. Bump
# them together at a release, with the version DOIs beside them, which no test
# can know before Zenodo issues them.
#
# Usage: cmake -DVERSION=<PROJECT_VERSION> -DSOURCE_DIR=<repository root>
#              -P check_versions.cmake

set(mismatches "")

# The one line of FILE matching REGEX, whose first group is the version.
function(check_version file regex)
    file(STRINGS "${SOURCE_DIR}/${file}" lines REGEX "${regex}")
    list(LENGTH lines count)
    if(NOT count EQUAL 1)
        set(mismatches "${mismatches}\n  ${file}: ${count} version lines, expected one"
            PARENT_SCOPE)
        return()
    endif()
    string(REGEX MATCH "${regex}" _ "${lines}")
    if(NOT CMAKE_MATCH_1 STREQUAL VERSION)
        set(mismatches "${mismatches}\n  ${file}: ${CMAKE_MATCH_1}" PARENT_SCOPE)
    endif()
endfunction()

check_version(pyproject.toml "^version = \"([^\"]*)\"$")
check_version(CITATION.cff "^version: \"([^\"]*)\"$")
check_version(README.md "^ +version += {([^}]*)},?$")
# The README's sentence naming the version DOI: "**v0.2.0 is [10.5281/...".
check_version(README.md "\\*\\*v([0-9.]+) is \\[")

if(mismatches)
    message(FATAL_ERROR "CMakeLists.txt's project(VERSION) is ${VERSION}, "
                        "but these say otherwise:${mismatches}")
endif()
