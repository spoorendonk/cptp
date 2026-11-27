# docs/options.md carries `cptp-solve --help` verbatim between two marker
# lines. MODE=check (the options_doc ctest test) fails when the block differs
# from what this build's cptp-solve prints; MODE=update (the
# cptp_options_doc_update target) rewrites the block in place, leaving the
# prose around it alone.
#
# Usage: cmake -DSOLVE=<cptp-solve> -DDOC=<docs/options.md> -DMODE=check|update
#              -P options_doc.cmake

set(begin_marker "<!-- cptp-solve --help: begin -->\n")
set(end_marker "<!-- cptp-solve --help: end -->")

execute_process(COMMAND "${SOLVE}" --help
    OUTPUT_VARIABLE help_stdout
    ERROR_VARIABLE help
    RESULT_VARIABLE help_result)
if(NOT help_result EQUAL 0)
    message(FATAL_ERROR "${SOLVE} --help exited with ${help_result}")
endif()
# The usage line names argv[0], which depends on how the binary was called.
string(REGEX REPLACE "Usage: [^\n]* <instance>" "Usage: cptp-solve <instance>"
       help "${help}")
set(block "${begin_marker}```text\n${help}```\n${end_marker}")

if(NOT EXISTS "${DOC}")
    message(FATAL_ERROR "${DOC} is missing")
endif()
file(READ "${DOC}" doc)
string(FIND "${doc}" "${begin_marker}" begin_pos)
string(FIND "${doc}" "${end_marker}" end_pos)
if(begin_pos EQUAL -1 OR end_pos EQUAL -1 OR end_pos LESS begin_pos)
    message(FATAL_ERROR "${DOC} lacks the marker lines around the help block")
endif()
string(LENGTH "${end_marker}" end_length)
math(EXPR block_length "${end_pos} + ${end_length} - ${begin_pos}")
string(SUBSTRING "${doc}" ${begin_pos} ${block_length} committed)

if(MODE STREQUAL "update")
    string(SUBSTRING "${doc}" 0 ${begin_pos} head)
    math(EXPR tail_pos "${begin_pos} + ${block_length}")
    string(SUBSTRING "${doc}" ${tail_pos} -1 tail)
    file(WRITE "${DOC}" "${head}${block}${tail}")
elseif(MODE STREQUAL "check")
    if(NOT committed STREQUAL block)
        # For the diff, beside the test (ctest runs it in the build tree).
        set(expected "${CMAKE_CURRENT_BINARY_DIR}/options_doc.expected")
        set(actual "${CMAKE_CURRENT_BINARY_DIR}/options_doc.committed")
        file(WRITE "${expected}" "${block}\n")
        file(WRITE "${actual}" "${committed}\n")
        execute_process(COMMAND diff -u "${actual}" "${expected}")
        message(FATAL_ERROR "${DOC} is out of date with cptp-solve --help: "
                            "build the cptp_options_doc_update target and commit it")
    endif()
else()
    message(FATAL_ERROR "MODE must be check or update, not '${MODE}'")
endif()
