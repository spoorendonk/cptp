if(NOT SKBUILD)
    # FindPython, not FindPython3: nanobind builds against `Python::Module`.
    find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
    include(FetchContent)
    FetchContent_Declare(nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        # 2.9+: stubgen types NumPy arrays as NDArray, not ArrayLike (#81).
        GIT_TAG v2.15.0
    )
    FetchContent_MakeAvailable(nanobind)
else()
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(nanobind CONFIG REQUIRED)
endif()

nanobind_add_module(cptp_python
    bindings.cpp
    NB_STATIC
)
set_target_properties(cptp_python PROPERTIES OUTPUT_NAME "_cptp")

target_link_libraries(cptp_python PRIVATE cptp)
target_compile_definitions(cptp_python PRIVATE CPTP_HAS_HIGHS)
target_include_directories(cptp_python PRIVATE
    ${highs_SOURCE_DIR}/src
    ${highs_BINARY_DIR}
)

# Only when the bindings were asked for (a wheel build sets this): built just
# for python_tests, they are not part of an install.
if(CPTP_BUILD_PYTHON)
    install(TARGETS cptp_python LIBRARY DESTINATION cptp)
endif()

# Type stubs (#81). The committed python/cptp/_cptp.pyi is what mypy and a wheel
# see; every build of the bindings regenerates it here, and the python_stub
# test fails when the two differ (refresh it with the cptp_python_stub_update
# target). A wheel build ships the committed copy and generates nothing.
#
# stubgen imports `cptp._cptp`, so the module is built into a `cptp` package
# whose `__init__.py` is empty: the real one imports numpy and the rest of the
# package, which stubgen needs none of.
if(NOT SKBUILD)
    set(CPTP_STUB_PACKAGE ${CMAKE_CURRENT_BINARY_DIR}/cptp)
    set(CPTP_STUB ${CPTP_STUB_PACKAGE}/_cptp.pyi)
    set(CPTP_STUB_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/cptp/_cptp.pyi)
    set_target_properties(cptp_python PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CPTP_STUB_PACKAGE})
    file(WRITE ${CPTP_STUB_PACKAGE}/__init__.py "")
    nanobind_add_stub(cptp_python_stub
        MODULE cptp._cptp
        OUTPUT ${CPTP_STUB}
        PYTHON_PATH ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS cptp_python)
    add_custom_target(cptp_python_stub_update
        COMMAND ${CMAKE_COMMAND} -E copy ${CPTP_STUB} ${CPTP_STUB_SOURCE}
        DEPENDS cptp_python_stub
        COMMENT "Refreshing python/cptp/_cptp.pyi from the build"
        VERBATIM)
    # For the python_stub test, which the top level registers.
    set(CPTP_STUB ${CPTP_STUB} PARENT_SCOPE)
    set(CPTP_STUB_SOURCE ${CPTP_STUB_SOURCE} PARENT_SCOPE)
endif()
