# Patch script for HiGHS: add HighsUserSeparator
# Called by FetchContent PATCH_COMMAND
# Idempotent: safe to run multiple times.

# Detect source layout: v1.13+ uses highs/ subdirectory
if(EXISTS "${SOURCE_DIR}/highs/mip")
    set(MIP_DIR "${SOURCE_DIR}/highs/mip")
    set(SRC_PREFIX "mip")
else()
    set(MIP_DIR "${SOURCE_DIR}/src/mip")
    set(SRC_PREFIX "mip")
endif()

# Copy HighsUserSeparator.h and .cpp
file(COPY "${PATCH_DIR}/HighsUserSeparator.h"
     DESTINATION "${MIP_DIR}/")
file(COPY "${PATCH_DIR}/HighsUserSeparator.cpp"
     DESTINATION "${MIP_DIR}/")

# ── 1. Patch cmake/sources.cmake: add HighsUserSeparator to source lists ──
file(READ "${SOURCE_DIR}/cmake/sources.cmake" SOURCES_CONTENT)

string(FIND "${SOURCES_CONTENT}" "HighsUserSeparator" _src_found)
if(_src_found EQUAL -1)
    # Add .cpp after HighsModkSeparator.cpp
    string(REPLACE
      "${SRC_PREFIX}/HighsModkSeparator.cpp"
      "${SRC_PREFIX}/HighsModkSeparator.cpp\n    ${SRC_PREFIX}/HighsUserSeparator.cpp"
      SOURCES_CONTENT "${SOURCES_CONTENT}")

    # Add .h after HighsModkSeparator.h
    string(REPLACE
      "${SRC_PREFIX}/HighsModkSeparator.h"
      "${SRC_PREFIX}/HighsModkSeparator.h\n    ${SRC_PREFIX}/HighsUserSeparator.h"
      SOURCES_CONTENT "${SOURCES_CONTENT}")

    file(WRITE "${SOURCE_DIR}/cmake/sources.cmake" "${SOURCES_CONTENT}")
    message(STATUS "Added HighsUserSeparator to cmake/sources.cmake")
else()
    message(STATUS "HighsUserSeparator already in cmake/sources.cmake, skipping")
endif()

# ── 2. Patch HighsSeparation.cpp: include + instantiate + lazy constraint at non-root ──
file(READ "${MIP_DIR}/HighsSeparation.cpp" CONTENT)

string(FIND "${CONTENT}" "HighsUserSeparator" _found)
if(_found EQUAL -1)
    # Add include for HighsUserSeparator.h after the last existing separator include
    string(REPLACE
      "#include \"mip/HighsTransformedLp.h\""
      "#include \"mip/HighsTransformedLp.h\"\n#include \"mip/HighsUserSeparator.h\""
      CONTENT "${CONTENT}")

    # Add HighsUserSeparator instantiation after the last separator emplace_back
    string(REPLACE
      "separators.emplace_back(new HighsModkSeparator(mipsolver));"
      "separators.emplace_back(new HighsModkSeparator(mipsolver));\n  separators.emplace_back(new HighsUserSeparator(mipsolver));"
      CONTENT "${CONTENT}")

    # Patch separate() entry condition: also enter separation when user callback
    # is registered and LP is integer feasible (lazy constraint at non-root nodes).
    # Original: if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())
    # New:      if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()))
    string(REPLACE
      "if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())"
      "if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()))"
      CONTENT "${CONTENT}")

    # Patch inner loop break: don't break on integer feasible if user callback is registered.
    # Original: if (ncuts == 0 || !lp->scaledOptimal(status) || lp->getFractionalIntegers().empty())
    # New:      if (ncuts == 0 || !lp->scaledOptimal(status) || (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback()))
    string(REPLACE
      "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          lp->getFractionalIntegers().empty())"
      "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback()))"
      CONTENT "${CONTENT}")

    file(WRITE "${MIP_DIR}/HighsSeparation.cpp" "${CONTENT}")
    message(STATUS "Applied HighsUserSeparator patch to HighsSeparation.cpp")
else()
    message(STATUS "HighsUserSeparator patch already applied, skipping")
endif()

# ── 3. Patch HighsMipSolverData.cpp: root separation + incumbent feasibility ──
file(READ "${MIP_DIR}/HighsMipSolverData.cpp" CONTENT2)

string(FIND "${CONTENT2}" "HighsUserSeparator" _found2)
if(_found2 EQUAL -1)
    # Add include for HighsUserSeparator.h
    string(REPLACE
      "#include \"mip/HighsMipSolverData.h\""
      "#include \"mip/HighsMipSolverData.h\"\n#include \"mip/HighsUserSeparator.h\""
      CONTENT2 "${CONTENT2}")

    # Modify the root separation loop condition to also run on integer solutions
    # when user separator has a callback (lazy constraint behavior).
    # Original: while (... && !lp.getFractionalIntegers().empty() && ...)
    # New: while (... && (!lp.getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()) && ...)
    string(REPLACE
      "!lp.getFractionalIntegers().empty() &&\n         stall < 3"
      "(!lp.getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()) &&\n         stall < 3"
      CONTENT2 "${CONTENT2}")

    # Patch addIncumbent(): reject solutions that fail user feasibility check.
    # Skip for sub-MIPs (different column space).
    # Original:  const bool execute_mip_solution_callback =
    # New:       if (!mipsolver.submip && !HighsUserSeparator::isFeasible(sol)) return false;
    #            const bool execute_mip_solution_callback =
    string(REPLACE
      "  const bool execute_mip_solution_callback ="
      "  // Reject solutions that violate user lazy constraints (e.g. SEC/GSEC)\n  if (!mipsolver.submip && !HighsUserSeparator::isFeasible(sol)) return false;\n  const bool execute_mip_solution_callback ="
      CONTENT2 "${CONTENT2}")

    # Patch runSetup(): also check user feasibility for MIP start solutions.
    # HiGHS accepts warm-start solutions directly without calling addIncumbent(),
    # bypassing our feasibility check. Add isFeasible() to the feasibility condition.
    # Original: if (feasible && solobj < upper_bound) {
    # New:      if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(incumbent);
    #           if (feasible && solobj < upper_bound) {
    string(REPLACE
      "    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;"
      "    if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(incumbent);\n    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;"
      CONTENT2 "${CONTENT2}")

    file(WRITE "${MIP_DIR}/HighsMipSolverData.cpp" "${CONTENT2}")
    message(STATUS "Applied lazy constraint + incumbent feasibility patch to HighsMipSolverData.cpp")
else()
    message(STATUS "HighsMipSolverData.cpp patch already applied, skipping")
endif()
