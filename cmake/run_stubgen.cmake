# nanobind's stubgen over this build's `cptp._cptp` (python/CMakeLists.txt),
# writing the type stubs the python_stub test compares with the committed
# python/cptp/_cptp.pyi.
#
# As in run_pytest.cmake, `python -S` with the site directories on PYTHONPATH,
# after the package the build made: no `.pth` file runs, so the editable install
# `uv sync` puts in .venv (whose import hook is consulted before PYTHONPATH, and
# serves the checkout's snapshot of the extension) cannot shadow the module the
# build just linked. The import is checked before the stub is written.
#
# Usage: cmake -DPYTHON=<exe> -DSTUBGEN=<stubgen.py> -DPACKAGE_DIR=<dir holding
#              cptp/> -DOUTPUT=<.pyi> -P run_stubgen.cmake

execute_process(
    COMMAND "${PYTHON}" -c
        "import site; print(*site.getsitepackages(), *([site.getusersitepackages()] if site.ENABLE_USER_SITE else []), sep='\\n')"
    OUTPUT_VARIABLE site_dirs
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${PYTHON} could not name its site-packages")
endif()
string(REPLACE "\n" ":" site_dirs "${site_dirs}")
set(env ${CMAKE_COMMAND} -E env "PYTHONPATH=${PACKAGE_DIR}:${site_dirs}")

execute_process(
    COMMAND ${env} "${PYTHON}" -S -c
        "import os, sys, cptp._cptp as m; d = os.path.realpath(os.path.dirname(m.__file__)); sys.exit(0 if d == os.path.realpath(sys.argv[1]) else 'cptp._cptp came from ' + d)"
        "${PACKAGE_DIR}/cptp"
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "stubgen would not read this build's cptp._cptp")
endif()

execute_process(
    COMMAND ${env} "${PYTHON}" -S "${STUBGEN}" -q -m cptp._cptp -o "${OUTPUT}"
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "stubgen failed (exit ${rc})")
endif()
