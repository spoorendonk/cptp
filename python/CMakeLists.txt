find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

include(FetchContent)
FetchContent_Declare(nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v2.4.0
)
FetchContent_MakeAvailable(nanobind)

nanobind_add_module(cptp_python
    bindings.cpp
    NB_STATIC
)
set_target_properties(cptp_python PROPERTIES OUTPUT_NAME "_cptp")
target_link_libraries(cptp_python PRIVATE cptp_model)
target_include_directories(cptp_python PRIVATE
    ${highs_SOURCE_DIR}/src
    ${highs_BINARY_DIR}
)

install(TARGETS cptp_python LIBRARY DESTINATION cptp)
