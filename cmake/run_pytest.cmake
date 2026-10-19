# The Python binding tests against this build's extension: the python_tests
# ctest test (CMakeLists.txt). The package is staged from the sources and the
# freshly built `_cptp` module, so the tests never see a stale extension.
#
# `python -S` with the interpreter's site-packages on PYTHONPATH: numpy and
# pytest import as usual, but no `.pth` file runs, so the editable install
# `uv sync` puts in .venv (whose import hook is consulted before PYTHONPATH,
# and serves the checkout's snapshot of the extension) cannot shadow the
# staged package.
#
# The CLI tests run this build's cptp-solve (CPTP_SOLVE), and the
# `python -m cptp` tests pass `-S` on to their child interpreter.
#
# Usage: cmake -DPYTHON=<exe> -DSOURCE_DIR=<root> -DSTAGE_DIR=<dir>
#              -DEXTENSION=<_cptp module> -DSOLVE=<cptp-solve> -P run_pytest.cmake

file(REMOVE_RECURSE "${STAGE_DIR}/cptp")
file(COPY "${SOURCE_DIR}/python/cptp" DESTINATION "${STAGE_DIR}"
    PATTERN "__pycache__" EXCLUDE)
file(COPY "${EXTENSION}" DESTINATION "${STAGE_DIR}/cptp")

execute_process(
    COMMAND "${PYTHON}" -c
        "import sysconfig; print(sysconfig.get_path('purelib')); print(sysconfig.get_path('platlib'))"
    OUTPUT_VARIABLE site_dirs
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${PYTHON} could not name its site-packages")
endif()
string(REPLACE "\n" ":" site_dirs "${site_dirs}")

execute_process(
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${STAGE_DIR}:${site_dirs}" "CPTP_SOLVE=${SOLVE}"
            "${PYTHON}" -S -m pytest -p no:cacheprovider "${SOURCE_DIR}/tests/python"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "pytest failed (exit ${rc})")
endif()
