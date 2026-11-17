# Patch script for HiGHS: add HighsUserSeparator + HighsUserPropagator
# Called by FetchContent PATCH_COMMAND (cmake/FetchHiGHS.cmake) with PATCH_HASH,
# a hash of this directory's patch files.
# Idempotent: a tree already patched with PATCH_HASH is left alone (the patch
# step runs on every configure); otherwise it restores the pristine HiGHS
# checkout, patches it and records PATCH_HASH.

if(NOT SOURCE_DIR OR NOT PATCH_DIR OR NOT PATCH_HASH OR NOT GIT_EXECUTABLE)
    message(FATAL_ERROR
      "apply_patch.cmake needs SOURCE_DIR, PATCH_DIR, PATCH_HASH and GIT_EXECUTABLE")
endif()

# Untracked, so the reset below removes it and a failed run leaves none.
set(_stamp "${SOURCE_DIR}/.cptp_patch_hash")
if(EXISTS "${_stamp}")
    file(READ "${_stamp}" _stamp_hash)
    if(_stamp_hash STREQUAL PATCH_HASH)
        message(STATUS "HiGHS already patched with ${PATCH_HASH}, skipping")
        return()
    endif()
endif()

# Start from the commit FetchContent checked out: drop an earlier run's edits
# (staged ones too) and copied sources, so every anchor below is matched
# against pristine HiGHS.
execute_process(COMMAND "${GIT_EXECUTABLE}" checkout HEAD -- .
    WORKING_DIRECTORY "${SOURCE_DIR}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${GIT_EXECUTABLE}" clean -fdq
    WORKING_DIRECTORY "${SOURCE_DIR}" COMMAND_ERROR_IS_FATAL ANY)

# Replace ANCHOR with REPLACEMENT in the variable VAR, which must hold exactly
# COUNT occurrences of ANCHOR. An anchor HiGHS no longer carries (a rename or a
# reformat after a bump) is a FATAL_ERROR here, not a silent no-op that shows
# up later at link time, or as a feature quietly off or a wrong optimum.
function(cptp_replace var count anchor replacement)
    string(LENGTH "${${var}}" _len)
    string(REPLACE "${anchor}" "" _stripped "${${var}}")
    string(LENGTH "${_stripped}" _len_stripped)
    string(LENGTH "${anchor}" _len_anchor)
    math(EXPR _found "(${_len} - ${_len_stripped}) / ${_len_anchor}")
    if(NOT _found EQUAL count)
        message(FATAL_ERROR
          "HiGHS patch anchor in ${var} matched ${_found} times, expected "
          "${count}. Update the anchor in apply_patch.cmake:\n${anchor}")
    endif()
    string(REPLACE "${anchor}" "${replacement}" _patched "${${var}}")
    set(${var} "${_patched}" PARENT_SCOPE)
endfunction()

# Source layout: v1.13+ keeps the sources under highs/
set(MIP_DIR "${SOURCE_DIR}/highs/mip")
set(SRC_PREFIX "mip")

# Copy HighsUserHooks.h, HighsUserSeparator and HighsUserPropagator (.h and
# .cpp). configure_file gives the copy a fresh mtime; file(COPY) would keep the
# source's, and a source older than its object (cp -p, rsync -t, tar) would not
# be rebuilt.
foreach(_name IN ITEMS HighsUserHooks.h
                       HighsUserSeparator.h HighsUserSeparator.cpp
                       HighsUserPropagator.h HighsUserPropagator.cpp)
    configure_file("${PATCH_DIR}/${_name}" "${MIP_DIR}/${_name}" COPYONLY)
endforeach()

# ── 1. Patch cmake/sources.cmake: add HighsUserSeparator to source lists ──
file(READ "${SOURCE_DIR}/cmake/sources.cmake" SOURCES_CONTENT)

# Add .cpp after HighsModkSeparator.cpp
cptp_replace(SOURCES_CONTENT 1
  "${SRC_PREFIX}/HighsModkSeparator.cpp"
  "${SRC_PREFIX}/HighsModkSeparator.cpp\n    ${SRC_PREFIX}/HighsUserSeparator.cpp\n    ${SRC_PREFIX}/HighsUserPropagator.cpp")

# Add .h after HighsModkSeparator.h
cptp_replace(SOURCES_CONTENT 1
  "${SRC_PREFIX}/HighsModkSeparator.h"
  "${SRC_PREFIX}/HighsModkSeparator.h\n    ${SRC_PREFIX}/HighsUserHooks.h\n    ${SRC_PREFIX}/HighsUserSeparator.h\n    ${SRC_PREFIX}/HighsUserPropagator.h")

# Keep every source this patch adds or edits out of HiGHS's unity build (on
# only if(APPLE) in highs/CMakeLists.txt): batched with neighbouring HiGHS
# sources, a file-scope name a hunk adds that a neighbour also has breaks only
# the macOS build, which no Linux gate builds. HiGHS does the same for
# io/filereaderlp/reader.cpp; sources.cmake is included from
# highs/CMakeLists.txt, which defines the targets, so the property applies.
# A hunk below that edits another HiGHS .cpp adds it to this list.
set(_cptp_patched_sources
  HighsUserSeparator.cpp HighsUserPropagator.cpp
  HighsSeparation.cpp HighsMipSolverData.cpp HighsSearch.cpp
  HighsMipSolver.cpp HighsMipWorker.cpp HighsPrimalHeuristics.cpp HighsCutPool.cpp)
list(TRANSFORM _cptp_patched_sources PREPEND "${SRC_PREFIX}/")
list(JOIN _cptp_patched_sources " " _cptp_patched_sources)
string(APPEND SOURCES_CONTENT
  "\n# cptp patch: keep the sources it adds or edits out of the unity build\n"
  "set_source_files_properties(${_cptp_patched_sources} PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)\n")

file(WRITE "${SOURCE_DIR}/cmake/sources.cmake" "${SOURCES_CONTENT}")
message(STATUS "Added HighsUserSeparator + HighsUserPropagator to cmake/sources.cmake")

# ── 2. Patch HighsSeparation.cpp: include + instantiate + lazy constraint at non-root ──
file(READ "${MIP_DIR}/HighsSeparation.cpp" CONTENT)

# Add include for HighsUserSeparator.h after the last existing separator include
cptp_replace(CONTENT 1
  "#include \"mip/HighsTransformedLp.h\""
  "#include \"mip/HighsTransformedLp.h\"\n#include \"mip/HighsUserSeparator.h\"")

# Add HighsUserSeparator instantiation after the last separator emplace_back
cptp_replace(CONTENT 1
  "separators.emplace_back(new HighsModkSeparator(mipsolver));"
  "separators.emplace_back(new HighsModkSeparator(mipsolver));\n  separators.emplace_back(new HighsUserSeparator(mipsolver));")

# Patch separate() entry condition: also enter separation when user callback
# is registered and LP is integer feasible (lazy constraint at non-root nodes).
# The hooks are the MIP solve's (HighsUserHooks), so a tree worker on another
# thread sees them too; `mipsolver` is lp->getMipSolver() here.
# Original: if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())
# New:      if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback(mipsolver)))
cptp_replace(CONTENT 1
  "if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())"
  "if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback(mipsolver)))")

# Patch inner loop break: don't break on integer feasible if user callback is registered.
# Original: if (ncuts == 0 || !lp->scaledOptimal(status) || lp->getFractionalIntegers().empty())
# New:      if (ncuts == 0 || !lp->scaledOptimal(status) || (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback(mipsolver)))
cptp_replace(CONTENT 1
  "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          lp->getFractionalIntegers().empty())"
  "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback(mipsolver)))")

file(WRITE "${MIP_DIR}/HighsSeparation.cpp" "${CONTENT}")
message(STATUS "Applied HighsUserSeparator patch to HighsSeparation.cpp")

# ── 3. Patch HighsMipSolverData.cpp: root separation + incumbent feasibility ──
file(READ "${MIP_DIR}/HighsMipSolverData.cpp" CONTENT2)

# Add include for HighsUserSeparator.h
cptp_replace(CONTENT2 1
  "#include \"mip/HighsMipSolverData.h\""
  "#include \"mip/HighsMipSolverData.h\"\n#include \"mip/HighsUserSeparator.h\"")

# Modify the root separation loop condition to also run on integer solutions
# when user separator has a callback (lazy constraint behavior).
# Original: while (... && !lp.getFractionalIntegers().empty() && ...)
# New: while (... && (!lp.getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()) && ...)
cptp_replace(CONTENT2 1
  "!getLp().getFractionalIntegers().empty() && stall < 3"
  "(!getLp().getFractionalIntegers().empty() || HighsUserSeparator::hasCallback(mipsolver)) && stall < 3")

# Patch addIncumbent(): reject solutions that fail user feasibility check.
# isFeasible(mipsolver, sol) expands a sub-MIP's reduced solution to its
# original space (postSolveStack.undoPrimal) before checking. Section 5 adds
# the same check to HighsMipWorker::addIncumbent, the path every incumbent of
# the parallel tree search takes while the parallel lock is active.
# Original:  const bool execute_mip_solution_callback =
# New:       if (!isFeasible(mipsolver, sol)) { mark rejected; return false; }
#            const bool execute_mip_solution_callback =
cptp_replace(CONTENT2 1
  "  const bool execute_mip_solution_callback ="
  "  // Reject solutions that violate user lazy constraints (e.g. SEC/GSEC)\n  if (!HighsUserSeparator::isFeasible(mipsolver, sol)) {\n    HighsUserSeparator::markIncumbentRejected();\n    return false;\n  }\n  const bool execute_mip_solution_callback =")

# Patch runSetup(): also check user feasibility for MIP start solutions.
# HiGHS accepts warm-start solutions directly without calling addIncumbent(),
# bypassing our feasibility check. Add isFeasible() to the feasibility condition.
# Original: if (feasible && solobj < upper_bound) {
# New:      if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(mipsolver, incumbent);
#           if (feasible && solobj < upper_bound) {
# The reassigned flag also gates the MIP-solution callback further down.
# From v1.16 checkAddSolution() declares it const: drop that so it can be
# reassigned.
cptp_replace(CONTENT2 1
  "    const bool feasible =\n        mipsolver.bound_violation_ <="
  "    bool feasible =\n        mipsolver.bound_violation_ <=")
cptp_replace(CONTENT2 1
  "    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;"
  "    if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(mipsolver, incumbent);\n    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;")

file(WRITE "${MIP_DIR}/HighsMipSolverData.cpp" "${CONTENT2}")
message(STATUS "Applied lazy constraint + incumbent feasibility patch to HighsMipSolverData.cpp")

# ── 4. Patch HighsSearch.cpp: add HighsUserPropagator call after reduced-cost fixing ──
file(READ "${MIP_DIR}/HighsSearch.cpp" CONTENT3)

# Add include for HighsUserPropagator.h
cptp_replace(CONTENT3 1
  "#include \"mip/HighsSearch.h\""
  "#include \"mip/HighsSearch.h\"\n#include \"mip/HighsUserPropagator.h\"")

# Insert propagator call after HighsRedcostFixing::propagateRedCost and
# before localdom.propagate(). The call is multi-line (parallel-MIP
# signature, since v1.15), so anchor on its closing arguments to stay unique.
cptp_replace(CONTENT3 1
  "mipworker.getPseudocost(), getUpperLimit());\n            localdom.propagate();"
  "mipworker.getPseudocost(), getUpperLimit());\n            HighsUserPropagator::propagate(localdom, mipsolver, *lp);\n            localdom.propagate();")

file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${CONTENT3}")
message(STATUS "Applied HighsUserPropagator patch to HighsSearch.cpp")

# ── 4c. Patch HighsSearch.cpp: keep integral nodes the lazy constraints reject ──
# evaluateNode() hands an integral LP solution to addIncumbent() and then prunes
# the node as bound exceeding. With lazy SECs the solution may contain subtours:
# section 3 rejects it, the upper limit does not move, and the pruned node can
# hold the optimum (M-n101-k10-97 proved -30918 instead of -32628 on v1.16).
# A rejected node stays open: branch() then branches on an unfixed integer
# column, and a node taken from the queue is separated (section 2).
file(READ "${MIP_DIR}/HighsSearch.cpp" CONTENT3C)

cptp_replace(CONTENT3C 1
  "#include \"mip/HighsSearch.h\""
  "#include \"mip/HighsSearch.h\"\n#include \"mip/HighsUserSeparator.h\"")

cptp_replace(CONTENT3C 1
  "          double cutoffbnd = getCutoffBound();\n          addIncumbent(lp->getLpSolver().getSolution().col_value,\n                       lp->getObjective(),\n                       inheuristic ? kSolutionSourceHeuristic\n                                   : kSolutionSourceEvaluateNode);\n          if (getUpperLimit() < cutoffbnd)\n            lp->setObjectiveLimit(getUpperLimit());\n\n          if (lp->unscaledDualFeasible(status)) {"
  "          double cutoffbnd = getCutoffBound();\n          HighsUserSeparator::takeIncumbentRejected();\n          addIncumbent(lp->getLpSolver().getSolution().col_value,\n                       lp->getObjective(),\n                       inheuristic ? kSolutionSourceHeuristic\n                                   : kSolutionSourceEvaluateNode);\n          const bool user_rejected =\n              HighsUserSeparator::takeIncumbentRejected();\n          if (getUpperLimit() < cutoffbnd)\n            lp->setObjectiveLimit(getUpperLimit());\n\n          if (!user_rejected && lp->unscaledDualFeasible(status)) {")

file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${CONTENT3C}")
message(STATUS "Applied lazy-rejection node patch to HighsSearch.cpp")

# ── 5. Hooks per MIP solve, so the parallel tree search's workers see them ──
# The setters fill a thread_local pending set on the Highs::run() caller.
# Tree workers run on HiGHS's executor threads (parallel/HighsTaskExecutor.h),
# so the callbacks live on the HighsMipSolver instead: a top-level run()
# captures the pending set before any task is spawned, and a sub-MIP (RINS,
# RENS), which may be started on any worker thread, copies its parent's.
file(READ "${MIP_DIR}/HighsMipSolver.h" CONTENT5H)
cptp_replace(CONTENT5H 1
  "struct HighsMipSolverData;"
  "struct HighsMipSolverData;\nstruct HighsUserHooks;")
cptp_replace(CONTENT5H 1
  "  void setParallelLock(bool lock) const;"
  "  void setParallelLock(bool lock) const;\n\n  // cptp patch: the user callbacks this solve runs (mip/HighsUserHooks.h)\n  std::shared_ptr<const HighsUserHooks> user_hooks_;")
file(WRITE "${MIP_DIR}/HighsMipSolver.h" "${CONTENT5H}")

file(READ "${MIP_DIR}/HighsMipSolver.cpp" CONTENT5)
cptp_replace(CONTENT5 1
  "#include \"mip/HighsMipSolver.h\""
  "#include \"mip/HighsMipSolver.h\"\n#include \"mip/HighsUserHooks.h\"")
# run(), not runMipPresolve(), which also creates a HighsMipSolverData.
cptp_replace(CONTENT5 1
  "fopen(options_mip_->mip_improving_solution_file.c_str(), \"w\");\n\n  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));"
  "fopen(options_mip_->mip_improving_solution_file.c_str(), \"w\");\n\n  if (!submip) user_hooks_ = HighsUserHooks::capture();\n  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));")
file(WRITE "${MIP_DIR}/HighsMipSolver.cpp" "${CONTENT5}")

file(READ "${MIP_DIR}/HighsPrimalHeuristics.cpp" CONTENT5P)
cptp_replace(CONTENT5P 1
  "  submipsolver.initialiseTerminator(mipsolver);"
  "  submipsolver.initialiseTerminator(mipsolver);\n  submipsolver.user_hooks_ = mipsolver.user_hooks_;")
file(WRITE "${MIP_DIR}/HighsPrimalHeuristics.cpp" "${CONTENT5P}")
message(STATUS "Applied per-solve user hooks patch to HighsMipSolver and HighsPrimalHeuristics")

# ── 6. Patch HighsMipWorker.cpp: lazy check on a tree worker's incumbent ──
# While the parallel lock is active, HighsSearch and HighsPrimalHeuristics hand
# incumbents to HighsMipWorker::addIncumbent, which tightens the worker's upper
# bound and limit and buffers the solution; HighsMipSolverData::addIncumbent
# (section 3) sees it only later, at syncSolutions. Check it here first, so a
# solution with subtours neither tightens the worker's limits nor lets
# evaluateNode prune its node (section 4c reads the rejection mark).
file(READ "${MIP_DIR}/HighsMipWorker.cpp" CONTENT6)
cptp_replace(CONTENT6 1
  "#include \"mip/HighsMipWorker.h\""
  "#include \"mip/HighsMipWorker.h\"\n#include \"mip/HighsUserSeparator.h\"")
cptp_replace(CONTENT6 1
  "                                  int solution_source) {\n  if (solobj < upper_bound) {"
  "                                  int solution_source) {\n  // Reject solutions that violate user lazy constraints (e.g. SEC/GSEC)\n  if (!HighsUserSeparator::isFeasible(mipsolver_, sol)) {\n    HighsUserSeparator::markIncumbentRejected();\n    return false;\n  }\n  if (solobj < upper_bound) {")
file(WRITE "${MIP_DIR}/HighsMipWorker.cpp" "${CONTENT6}")
message(STATUS "Applied worker incumbent feasibility patch to HighsMipWorker.cpp")

# ── 7. Patch HighsCutPool.cpp: keep a cut's propagate flag through a sync ──
# syncCutPool copies a worker's cuts into the global pool with addCut's default
# propagate = true, so a cut added with propagate = false (cptp's RCI cuts)
# would be propagated once synced. columnsLinked(i) is the cut's current flag:
# false also for a cut addCut found too dense, or one the worker pool later
# unlinked to stay within its propagation budget. Such a cut syncs as
# non-propagated, which is weaker but never wrong: it stays a valid cut.
file(READ "${MIP_DIR}/HighsCutPool.cpp" CONTENT7)
cptp_replace(CONTENT7 1
  "      syncpool.addCut(mipsolver, idxs.data(), vals.data(), Rlen, rhs_[i],\n                      rowintegral[i]);"
  "      syncpool.addCut(mipsolver, idxs.data(), vals.data(), Rlen, rhs_[i],\n                      rowintegral[i], matrix_.columnsLinked(i));")
file(WRITE "${MIP_DIR}/HighsCutPool.cpp" "${CONTENT7}")
message(STATUS "Applied cut propagate-flag sync patch to HighsCutPool.cpp")

file(WRITE "${_stamp}" "${PATCH_HASH}")
