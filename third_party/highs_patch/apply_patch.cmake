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

# Copy HighsUserSeparator and HighsUserPropagator (.h and .cpp). configure_file
# gives the copy a fresh mtime; file(COPY) would keep the source's, and a
# source older than its object (cp -p, rsync -t, tar) would not be rebuilt.
foreach(_name IN ITEMS HighsUserSeparator.h HighsUserSeparator.cpp
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
  "${SRC_PREFIX}/HighsModkSeparator.h\n    ${SRC_PREFIX}/HighsUserSeparator.h\n    ${SRC_PREFIX}/HighsUserPropagator.h")

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
  HighsLpRelaxation.cpp)
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
# Original: if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())
# New:      if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()))
cptp_replace(CONTENT 1
  "if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty())"
  "if (lp->scaledOptimal(status) && (!lp->getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()))")

# Patch inner loop break: don't break on integer feasible if user callback is registered.
# Original: if (ncuts == 0 || !lp->scaledOptimal(status) || lp->getFractionalIntegers().empty())
# New:      if (ncuts == 0 || !lp->scaledOptimal(status) || (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback()))
cptp_replace(CONTENT 1
  "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          lp->getFractionalIntegers().empty())"
  "if (ncuts == 0 || !lp->scaledOptimal(status) ||\n          (lp->getFractionalIntegers().empty() && !HighsUserSeparator::hasCallback()))")

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
  "(!getLp().getFractionalIntegers().empty() || HighsUserSeparator::hasCallback()) && stall < 3")

# Patch addIncumbent(): reject solutions that fail user feasibility check.
# For sub-MIPs (with presolve on), expand reduced solution to original space
# via postSolveStack.undoPrimal before checking feasibility.
# Original:  const bool execute_mip_solution_callback =
# New:       { bool user_feasible = true;
#              if (mipsolver.submip) { expand via undoPrimal, check }
#              else { check directly }
#              if (!user_feasible) { mark rejected; return false; } }
#            const bool execute_mip_solution_callback =
cptp_replace(CONTENT2 1
  "  const bool execute_mip_solution_callback ="
  "  // Reject solutions that violate user lazy constraints (e.g. SEC/GSEC)\n  {\n    bool user_feasible;\n    if (mipsolver.submip) {\n      // Sub-MIP: sol is in reduced (presolved) space; expand to original\n      HighsSolution temp_sol;\n      temp_sol.col_value = sol;\n      temp_sol.value_valid = true;\n      postSolveStack.undoPrimal(*mipsolver.options_mip_, temp_sol);\n      user_feasible = HighsUserSeparator::isFeasible(temp_sol.col_value);\n    } else {\n      user_feasible = HighsUserSeparator::isFeasible(sol);\n    }\n    if (!user_feasible) {\n      HighsUserSeparator::markIncumbentRejected();\n      return false;\n    }\n  }\n  const bool execute_mip_solution_callback =")

# Patch runSetup(): also check user feasibility for MIP start solutions.
# HiGHS accepts warm-start solutions directly without calling addIncumbent(),
# bypassing our feasibility check. Add isFeasible() to the feasibility condition.
# Original: if (feasible && solobj < upper_bound) {
# New:      if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(incumbent);
#           if (feasible && solobj < upper_bound) {
# The reassigned flag also gates the MIP-solution callback further down.
# From v1.16 checkAddSolution() declares it const: drop that so it can be
# reassigned.
cptp_replace(CONTENT2 1
  "    const bool feasible =\n        mipsolver.bound_violation_ <="
  "    bool feasible =\n        mipsolver.bound_violation_ <=")
cptp_replace(CONTENT2 1
  "    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;"
  "    if (feasible && !mipsolver.submip) feasible = HighsUserSeparator::isFeasible(incumbent);\n    if (feasible && solobj < upper_bound) {\n      double prev_upper_bound = upper_bound;")

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

# ── 4a. Patch HighsLpRelaxation.h: add kBranching origin + methods ──
file(READ "${MIP_DIR}/HighsLpRelaxation.h" LPRELAX_H)

# Add kBranching to the Origin enum (match the full enum context to avoid
# accidentally replacing kCutPool in factory methods)
cptp_replace(LPRELAX_H 1
  "      kCutPool,\n    };"
  "      kCutPool,\n      kBranching,\n    };")

# Add static factory method branching() after model() factory
cptp_replace(LPRELAX_H 1
  "static LpRow model(HighsInt index) { return LpRow{kModel, index, 0, -1}; }"
  "static LpRow model(HighsInt index) { return LpRow{kModel, index, 0, -1}; }\n    static LpRow branching() { return LpRow{kBranching, -1, 0, -1}; }")

# Add method declarations after removeObsoleteRows (unique, appears once)
cptp_replace(LPRELAX_H 1
  "void removeObsoleteRows(bool notifyPool = true);"
  "void removeObsoleteRows(bool notifyPool = true);\n\n  HighsInt addBranchingRow(HighsInt nnz, const HighsInt* inds,\n                            const double* vals, double lo, double hi);\n\n  void removeLastBranchingRow();")

file(WRITE "${MIP_DIR}/HighsLpRelaxation.h" "${LPRELAX_H}")
message(STATUS "Applied kBranching patch to HighsLpRelaxation.h")

# ── 4b. Patch HighsLpRelaxation.cpp: implement addBranchingRow/removeLastBranchingRow ──
file(READ "${MIP_DIR}/HighsLpRelaxation.cpp" LPRELAX_CPP)

# Handle kBranching in get(): return len=0 (branching rows have no cut/model data)
cptp_replace(LPRELAX_CPP 1
  "  switch (origin) {\n    case kCutPool:\n      mipsolver.mipdata_->cutpools[cutpoolindex].getCut(index, len, inds, vals);\n      break;\n    case kModel:\n      mipsolver.mipdata_->getRow(index, len, inds, vals);\n  };"
  "  switch (origin) {\n    case kCutPool:\n      mipsolver.mipdata_->cutpools[cutpoolindex].getCut(index, len, inds, vals);\n      break;\n    case kModel:\n      mipsolver.mipdata_->getRow(index, len, inds, vals);\n      break;\n    case kBranching:\n      len = 0; inds = nullptr; vals = nullptr;\n      break;\n  };")

# Handle kBranching in getRowLen(): return 0
cptp_replace(LPRELAX_CPP 1
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].getRowLength(index);\n    case kModel:\n      return mipsolver.mipdata_->ARstart_[index + 1] -\n             mipsolver.mipdata_->ARstart_[index];\n  };"
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].getRowLength(index);\n    case kModel:\n      return mipsolver.mipdata_->ARstart_[index + 1] -\n             mipsolver.mipdata_->ARstart_[index];\n    case kBranching:\n      return 0;\n  };")

# Handle kBranching in isIntegral(): return false
cptp_replace(LPRELAX_CPP 1
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].cutIsIntegral(index);\n    case kModel:\n      return mipsolver.mipdata_->rowintegral[index];\n  };"
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].cutIsIntegral(index);\n    case kModel:\n      return mipsolver.mipdata_->rowintegral[index];\n    case kBranching:\n      return false;\n  };")

# Handle kBranching in getMaxAbsVal(): return 0
cptp_replace(LPRELAX_CPP 1
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].getMaxAbsCutCoef(index);\n    case kModel:\n      return mipsolver.mipdata_->maxAbsRowCoef[index];\n  };"
  "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpools[cutpoolindex].getMaxAbsCutCoef(index);\n    case kModel:\n      return mipsolver.mipdata_->maxAbsRowCoef[index];\n    case kBranching:\n      return 0.0;\n  };")

# Handle kBranching in slackLower(): return row lower bound
cptp_replace(LPRELAX_CPP 1
  "  switch (lprows[row].origin) {\n    case LpRow::kCutPool:\n      return globaldom.getMinCutActivity(\n          mipsolver.mipdata_->cutpools[lprows[row].cutpoolindex],\n          lprows[row].index);\n    case LpRow::kModel:"
  "  switch (lprows[row].origin) {\n    case LpRow::kBranching:\n      return rowLower(row);\n    case LpRow::kCutPool:\n      return globaldom.getMinCutActivity(\n          mipsolver.mipdata_->cutpools[lprows[row].cutpoolindex],\n          lprows[row].index);\n    case LpRow::kModel:")

# Handle kBranching in slackUpper(): return rowupper
cptp_replace(LPRELAX_CPP 1
  "  switch (lprows[row].origin) {\n    case LpRow::kCutPool:\n      return rowupper;\n    case LpRow::kModel:"
  "  switch (lprows[row].origin) {\n    case LpRow::kBranching:\n      return rowupper;\n    case LpRow::kCutPool:\n      return rowupper;\n    case LpRow::kModel:")

# Append addBranchingRow and removeLastBranchingRow implementations at the end
string(APPEND LPRELAX_CPP "
// ── Dynamic branching rows ──

HighsInt HighsLpRelaxation::addBranchingRow(
    HighsInt nnz, const HighsInt* inds, const double* vals,
    double lo, double hi) {
  status = Status::kNotSet;
  currentbasisstored = false;
  basischeckpoint.reset();
  HighsInt row = lprows.size();
  lprows.push_back(LpRow::branching());
  lpsolver.addRow(lo, hi, nnz, inds, vals);
  return row;
}

void HighsLpRelaxation::removeLastBranchingRow() {
  for (HighsInt i = (HighsInt)lprows.size() - 1; i >= 0; --i) {
    if (lprows[i].origin == LpRow::Origin::kBranching) {
      lpsolver.deleteRows(i, i);
      lprows.erase(lprows.begin() + i);
      status = Status::kNotSet;
      currentbasisstored = false;
      basischeckpoint.reset();
      return;
    }
  }
}
")

file(WRITE "${MIP_DIR}/HighsLpRelaxation.cpp" "${LPRELAX_CPP}")
message(STATUS "Applied addBranchingRow/removeLastBranchingRow patch to HighsLpRelaxation.cpp")

# ── 5a. Patch HighsSearch.h: add hp_stack_pos and hp_data_idx to NodeData ──
file(READ "${MIP_DIR}/HighsSearch.h" SEARCH_H)

# Add fields to NodeData struct, after domgchgStackPos
cptp_replace(SEARCH_H 1
  "HighsInt domgchgStackPos;"
  "HighsInt domgchgStackPos;\n    HighsInt hp_stack_pos = 0;\n    HighsInt hp_data_idx = -1;")

file(WRITE "${MIP_DIR}/HighsSearch.h" "${SEARCH_H}")
message(STATUS "Applied hp_stack_pos/hp_data_idx patch to HighsSearch.h")

# ── 5b. Patch HighsSearch.cpp: hyperplane branching in branch() ──
file(READ "${MIP_DIR}/HighsSearch.cpp" SEARCH_CONTENT)

# Add includes
cptp_replace(SEARCH_CONTENT 1
  "#include \"mip/HighsSearch.h\""
  "#include \"mip/HighsSearch.h\"\n#include \"mip/HighsUserSeparator.h\"")

# Insert hyperplane candidate evaluation + dynamic row branching
# before normal child creation.
cptp_replace(SEARCH_CONTENT 1
  "  // finally open a new node with the branching decision added\n  // and remember that we have one open subtree left\n  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  currnode.opensubtrees = 1;\n  nodestack.emplace_back(\n      std::max(childLb, currnode.lower_bound), currnode.estimate,\n      currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n  nodestack.back().domgchgStackPos = domchgPos;\n\n  return NodeResult::kBranched;"
  "  // Hyperplane branching: evaluate candidates vs variable branching\n  const auto& branchCb = HighsUserSeparator::getBranchingCallback();\n  if (branchCb && !mipsolver.submip) {\n    auto candidates = branchCb(mipsolver);\n    if (!candidates.empty()) {\n      // Use pseudocost estimates for var_score (comparable with HP scoring).\n      // SB-based childLb/other_child_lb are only meaningful when minreliable > 0.\n      double var_score = 0.0;\n      if (currnode.branchingdecision.column >= 0) {\n        double vd = pseudocost.getPseudocostDown(\n            currnode.branchingdecision.column, currnode.branching_point);\n        double vu = pseudocost.getPseudocostUp(\n            currnode.branchingdecision.column, currnode.branching_point);\n        var_score = std::min(vd, vu) * std::max(vd, vu);\n      }\n\n      int winner = HighsUserSeparator::evaluateCandidates(\n          candidates, var_score, mipsolver, (int)getCurrentDepth());\n      if (winner >= 0) {\n        auto& cand = candidates[winner];\n        const auto& sol = lp->getSolution().col_value;\n        double lhs = 0.0;\n        for (int j = 0; j < (int)cand.indices.size(); ++j)\n          lhs += cand.values[j] * sol[cand.indices[j]];\n\n        // Store hyperplane data for flip\n        HighsInt data_idx = HighsUserSeparator::allocHyperplaneData();\n        auto& data = HighsUserSeparator::getHyperplaneData(data_idx);\n        data.indices = cand.indices;\n        data.values = cand.values;\n        data.branch_value = lhs;\n        currnode.hp_data_idx = data_idx;\n\n        // Add branching row (down branch: lhs <= floor)\n        HighsInt hpPos = HighsUserSeparator::stackSize();\n        HighsUserSeparator::pushBranching(\n            mipsolver.mipdata_->getLp(), data,\n            -kHighsInf, std::floor(lhs));\n\n        // Push no-op domain marker so backtrack() has a kBranching sentinel\n        HighsInt domchgPos = localdom.getDomainChangeStack().size();\n        localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n        currnode.opensubtrees = 1;\n        nodestack.emplace_back(\n            currnode.lower_bound, currnode.estimate,\n            currnode.nodeBasis, nullptr);\n        nodestack.back().domgchgStackPos = domchgPos;\n        nodestack.back().hp_stack_pos = hpPos;\n\n        return NodeResult::kBranched;\n      }\n    }\n  }\n\n  // finally open a new node with the branching decision added\n  // and remember that we have one open subtree left\n  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  currnode.opensubtrees = 1;\n  nodestack.emplace_back(\n      std::max(childLb, currnode.lower_bound), currnode.estimate,\n      currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n  nodestack.back().domgchgStackPos = domchgPos;\n  nodestack.back().hp_stack_pos = HighsUserSeparator::stackSize();\n\n  return NodeResult::kBranched;")

# ── 5c. Patch backtrack pop_back + localdom.backtrack(): restore branching rows ──
# The pattern has #ifndef NDEBUG around the backtrack() return value: match
# nodestack.pop_back() followed by that #ifndef block. It appears twice
# (backtrack, backtrackPlunge); backtrackUntilDepth has its own below.
cptp_replace(SEARCH_CONTENT 2
  "      nodestack.pop_back();\n#ifndef NDEBUG\n      HighsDomainChange branchchg =\n#endif\n          localdom.backtrack();"
  "      HighsInt hpPos_pop = nodestack.back().hp_stack_pos;\n      nodestack.pop_back();\n#ifndef NDEBUG\n      HighsDomainChange branchchg =\n#endif\n          localdom.backtrack();\n      HighsUserSeparator::restoreBranching(mipsolver.mipdata_->getLp(), hpPos_pop);")

# backtrackUntilDepth has a slightly different pattern (4-space indent)
cptp_replace(SEARCH_CONTENT 1
  "    nodestack.pop_back();\n\n#ifndef NDEBUG\n    HighsDomainChange branchchg =\n#endif\n        localdom.backtrack();"
  "    HighsInt hpPos_pop = nodestack.back().hp_stack_pos;\n    nodestack.pop_back();\n\n#ifndef NDEBUG\n    HighsDomainChange branchchg =\n#endif\n        localdom.backtrack();\n    HighsUserSeparator::restoreBranching(mipsolver.mipdata_->getLp(), hpPos_pop);")

# ── 5d. Patch flip sections in backtrack() and backtrackPlunge() ──
# When hp_data_idx >= 0, we skip the branchingdecision flip, push a no-op
# marker, and add the up-branch hyperplane row.
# This pattern appears twice (backtrack + backtrackPlunge) with 4-space indent;
# from v1.16 both declare numChangedCols as size_t.
cptp_replace(SEARCH_CONTENT 2
  "    size_t numChangedCols = localdom.getChangedCols().size();\n    bool passStabilizerToChildNode =\n        orbitsValidInChildNode(currnode.branchingdecision);\n    localdom.changeBound(currnode.branchingdecision);"
  "    size_t numChangedCols = localdom.getChangedCols().size();\n    HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n    if (currnode.hp_data_idx >= 0) {\n      // Hyperplane flip: push no-op marker + up-branch row\n      auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n      localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n      HighsUserSeparator::pushBranching(\n          mipsolver.mipdata_->getLp(), flipData,\n          std::ceil(flipData.branch_value), kHighsInf);\n    } else {\n      localdom.changeBound(currnode.branchingdecision);\n    }\n    bool passStabilizerToChildNode =\n        (currnode.hp_data_idx < 0) && orbitsValidInChildNode(currnode.branchingdecision);")

# ── 5e. Patch flip in backtrackUntilDepth() ──
# This one has 2-space indent.
cptp_replace(SEARCH_CONTENT 1
  "  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  nodestack.emplace_back(\n      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n\n  lp->flushDomain(localdom);\n  nodestack.back().domgchgStackPos = domchgPos;"
  "  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n  HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n  if (currnode.hp_data_idx >= 0) {\n    auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n    localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n    HighsUserSeparator::pushBranching(\n        mipsolver.mipdata_->getLp(), flipData,\n        std::ceil(flipData.branch_value), kHighsInf);\n  } else {\n    localdom.changeBound(currnode.branchingdecision);\n  }\n  bool passStabilizerToChildNode =\n      (currnode.hp_data_idx < 0) && orbitsValidInChildNode(currnode.branchingdecision);\n  nodestack.emplace_back(\n      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n\n  lp->flushDomain(localdom);\n  nodestack.back().domgchgStackPos = domchgPos;\n  nodestack.back().hp_stack_pos = hpPos_flip;")

# ── 5f. Set hp_stack_pos on flip child in backtrack() and backtrackPlunge() ──
# After the flip's emplace_back, set hp_stack_pos. The pattern ends with:
#   lp->flushDomain(localdom);
#   nodestack.back().domgchgStackPos = domchgPos;
#   break;
# This appears in both backtrack() and backtrackPlunge().
cptp_replace(SEARCH_CONTENT 2
  "    lp->flushDomain(localdom);\n    nodestack.back().domgchgStackPos = domchgPos;\n    break;"
  "    lp->flushDomain(localdom);\n    nodestack.back().domgchgStackPos = domchgPos;\n    nodestack.back().hp_stack_pos = hpPos_flip;\n    break;")

# ── 5g. Restore branching rows when flip child is pruned or queued ──
# When the flip child is pruned (infeasible/cutoff) or sent to the node queue,
# localdom.backtrack() is called but branching rows must also be removed.
# The prune pattern is: "localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight)"
# This appears in backtrack() and backtrackPlunge().
cptp_replace(SEARCH_CONTENT 2
  "      localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());\n      continue;"
  "      HighsUserSeparator::restoreBranching(mipsolver.mipdata_->getLp(), hpPos_flip);\n      localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());\n      continue;")

# ── 5h. Prevent hyperplane-branched nodes from going to the node queue ──
# The node queue stores only domain changes, not LP branching rows. The guard
# sits right before the nodeToQueue block (v1.16 has no ancestor score loop
# between that and the flag's definition any more).
cptp_replace(SEARCH_CONTENT 1
  "    if (nodeToQueue) {\n      // if (!mipsolver.submip) printf(\"node goes to queue\\n\");"
  "    if (currnode.hp_data_idx >= 0) nodeToQueue = false;  // hyperplane nodes can't be queued\n    if (nodeToQueue) {\n      // if (!mipsolver.submip) printf(\"node goes to queue\\n\");")

# ── 5i. Suppress debug asserts for hyperplane branches ──
# The asserts check branchchg matches branchingdecision. For hyperplane branches
# the domain change is a no-op marker, not the branchingdecision.
# Guard with hp_data_idx check.

# Pattern in backtrack() and backtrackPlunge() (6-space indent, multi-line assert)
cptp_replace(SEARCH_CONTENT 2
  "      assert(\n          (branchchg.boundtype == HighsBoundType::kLower &&\n           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n          (branchchg.boundtype == HighsBoundType::kUpper &&\n           branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n      assert(branchchg.boundtype ==\n             nodestack.back().branchingdecision.boundtype);\n      assert(branchchg.column == nodestack.back().branchingdecision.column);"
  "      if (nodestack.back().hp_data_idx < 0) {\n      assert(\n          (branchchg.boundtype == HighsBoundType::kLower &&\n           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n          (branchchg.boundtype == HighsBoundType::kUpper &&\n           branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n      assert(branchchg.boundtype ==\n             nodestack.back().branchingdecision.boundtype);\n      assert(branchchg.column == nodestack.back().branchingdecision.column);\n      }")

# Pattern in backtrackUntilDepth() (4-space indent, single-line asserts)
cptp_replace(SEARCH_CONTENT 1
  "    assert(\n        (branchchg.boundtype == HighsBoundType::kLower &&\n         branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n        (branchchg.boundtype == HighsBoundType::kUpper &&\n         branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);\n    assert(branchchg.column == nodestack.back().branchingdecision.column);"
  "    if (nodestack.back().hp_data_idx < 0) {\n    assert(\n        (branchchg.boundtype == HighsBoundType::kLower &&\n         branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n        (branchchg.boundtype == HighsBoundType::kUpper &&\n         branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);\n    assert(branchchg.column == nodestack.back().branchingdecision.column);\n    }")

file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${SEARCH_CONTENT}")
message(STATUS "Applied hyperplane branching patch to HighsSearch.cpp")

file(WRITE "${_stamp}" "${PATCH_HASH}")
