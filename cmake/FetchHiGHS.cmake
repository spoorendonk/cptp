include(FetchContent)

# Disable HiGHS components we don't need
set(HIGHS_NO_DEFAULT_THREADS ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PDLP OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

find_package(Git REQUIRED)

# FetchContent re-runs PATCH_COMMAND only when the command line changes, so it
# carries a hash of every file in the patch directory: a change to the script
# or a patch source re-runs it on the next configure, and apply_patch.cmake
# starts from the pristine checkout. CONFIGURE_DEPENDS on the glob and each
# file makes a plain build re-run configure after such a change.
set(_cptp_highs_patch_dir ${CMAKE_CURRENT_SOURCE_DIR}/third_party/highs_patch)
file(GLOB _cptp_highs_patch_files LIST_DIRECTORIES false CONFIGURE_DEPENDS
    ${_cptp_highs_patch_dir}/*)
list(SORT _cptp_highs_patch_files)
set(_cptp_highs_patch_digests "")
foreach(_file IN LISTS _cptp_highs_patch_files)
    file(SHA256 ${_file} _digest)
    cmake_path(GET _file FILENAME _name)
    string(APPEND _cptp_highs_patch_digests "${_name} ${_digest}\n")
endforeach()
string(SHA256 _cptp_highs_patch_hash "${_cptp_highs_patch_digests}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_cptp_highs_patch_files})

FetchContent_Declare(highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG        v1.16.0
    PATCH_COMMAND ${CMAKE_COMMAND}
        -DPATCH_DIR=${_cptp_highs_patch_dir}
        -DPATCH_HASH=${_cptp_highs_patch_hash}
        -DSOURCE_DIR=<SOURCE_DIR>
        -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
        -P ${_cptp_highs_patch_dir}/apply_patch.cmake
)

FetchContent_MakeAvailable(highs)
