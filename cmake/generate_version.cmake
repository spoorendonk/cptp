# Writes version.h with the current git describe string.
#
# Called from a custom target at build time (not configure time), so the hash
# stays in sync with HEAD without re-running cmake. configure_file only touches
# the output when content changes, so downstream rebuilds are minimized.
#
# Expected -D variables: SRC_DIR, CPTP_VERSION, VERSION_IN, VERSION_OUT.

find_package(Git QUIET)
set(CPTP_GIT_DESCRIBE "unknown")
if(GIT_FOUND)
    execute_process(
        # --long forces the -<n>-g<sha> suffix even on an exact tag, so the
        # commit sha is always present in --version output (not just for
        # off-tag/dirty builds). In a tagless tree, --always falls back to a
        # bare sha and the long form is dropped.
        COMMAND ${GIT_EXECUTABLE} describe --tags --long --always --dirty --abbrev=7
        WORKING_DIRECTORY ${SRC_DIR}
        OUTPUT_VARIABLE _describe
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _rc)
    if(_rc EQUAL 0 AND _describe)
        set(CPTP_GIT_DESCRIBE "${_describe}")
    endif()
endif()
configure_file(${VERSION_IN} ${VERSION_OUT} @ONLY)
