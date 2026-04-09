include(FetchContent)

# Disable HiGHS components we don't need
set(HIGHS_NO_DEFAULT_THREADS ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PDLP OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

FetchContent_Declare(highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG        v1.14.0
    PATCH_COMMAND ${CMAKE_COMMAND}
        -DPATCH_DIR=${CMAKE_CURRENT_SOURCE_DIR}/third_party/highs_patch
        -DSOURCE_DIR=<SOURCE_DIR>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/third_party/highs_patch/apply_patch.cmake
)

FetchContent_MakeAvailable(highs)
