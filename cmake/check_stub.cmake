# The python_stub ctest test (CMakeLists.txt): the committed type stubs for
# cptp._cptp match those nanobind's stubgen wrote for this build's extension
# (python/CMakeLists.txt), so mypy checks the package against what the
# bindings really export.
#
# Usage: cmake -DGENERATED=<built .pyi> -DCOMMITTED=<python/cptp/_cptp.pyi>
#              -P check_stub.cmake

if(NOT EXISTS "${GENERATED}")
    message(FATAL_ERROR "${GENERATED} was not generated: build the "
                        "cptp_python_stub target")
endif()
if(NOT EXISTS "${COMMITTED}")
    message(FATAL_ERROR "${COMMITTED} is missing: build the "
                        "cptp_python_stub_update target and commit it")
endif()
file(READ "${GENERATED}" generated)
file(READ "${COMMITTED}" committed)
if(NOT generated STREQUAL committed)
    execute_process(COMMAND diff -u "${COMMITTED}" "${GENERATED}")
    message(FATAL_ERROR "${COMMITTED} is out of date with the bindings: "
                        "build the cptp_python_stub_update target and commit it")
endif()
