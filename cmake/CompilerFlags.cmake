# Compiler flags applied only to our own targets (not dependencies).
# Use target_compile_options in CMakeLists.txt instead of global flags.

function(cptp_set_compiler_flags target)
    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic
        -Wno-unused-parameter
        $<$<CONFIG:Release>:-O3 -march=native>
    )
endfunction()
