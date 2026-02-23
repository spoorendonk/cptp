# Patch script for HiGHS: add HighsUserSeparator + HighsUserPropagator
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

# Copy HighsUserPropagator.h and .cpp
file(COPY "${PATCH_DIR}/HighsUserPropagator.h"
     DESTINATION "${MIP_DIR}/")
file(COPY "${PATCH_DIR}/HighsUserPropagator.cpp"
     DESTINATION "${MIP_DIR}/")

# ── 1. Patch cmake/sources.cmake: add HighsUserSeparator to source lists ──
file(READ "${SOURCE_DIR}/cmake/sources.cmake" SOURCES_CONTENT)

string(FIND "${SOURCES_CONTENT}" "HighsUserSeparator" _src_found)
if(_src_found EQUAL -1)
    # Add .cpp after HighsModkSeparator.cpp
    string(REPLACE
      "${SRC_PREFIX}/HighsModkSeparator.cpp"
      "${SRC_PREFIX}/HighsModkSeparator.cpp\n    ${SRC_PREFIX}/HighsUserSeparator.cpp\n    ${SRC_PREFIX}/HighsUserPropagator.cpp"
      SOURCES_CONTENT "${SOURCES_CONTENT}")

    # Add .h after HighsModkSeparator.h
    string(REPLACE
      "${SRC_PREFIX}/HighsModkSeparator.h"
      "${SRC_PREFIX}/HighsModkSeparator.h\n    ${SRC_PREFIX}/HighsUserSeparator.h\n    ${SRC_PREFIX}/HighsUserPropagator.h"
      SOURCES_CONTENT "${SOURCES_CONTENT}")

    file(WRITE "${SOURCE_DIR}/cmake/sources.cmake" "${SOURCES_CONTENT}")
    message(STATUS "Added HighsUserSeparator + HighsUserPropagator to cmake/sources.cmake")
else()
    # HighsUserSeparator present; check if HighsUserPropagator also present
    string(FIND "${SOURCES_CONTENT}" "HighsUserPropagator" _prop_found)
    if(_prop_found EQUAL -1)
        string(REPLACE
          "${SRC_PREFIX}/HighsUserSeparator.cpp"
          "${SRC_PREFIX}/HighsUserSeparator.cpp\n    ${SRC_PREFIX}/HighsUserPropagator.cpp"
          SOURCES_CONTENT "${SOURCES_CONTENT}")
        string(REPLACE
          "${SRC_PREFIX}/HighsUserSeparator.h"
          "${SRC_PREFIX}/HighsUserSeparator.h\n    ${SRC_PREFIX}/HighsUserPropagator.h"
          SOURCES_CONTENT "${SOURCES_CONTENT}")
        file(WRITE "${SOURCE_DIR}/cmake/sources.cmake" "${SOURCES_CONTENT}")
        message(STATUS "Added HighsUserPropagator to cmake/sources.cmake")
    else()
        message(STATUS "HighsUserSeparator + HighsUserPropagator already in cmake/sources.cmake, skipping")
    endif()
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
    # For sub-MIPs (with presolve on), expand reduced solution to original space
    # via postSolveStack.undoPrimal before checking feasibility.
    # Original:  const bool execute_mip_solution_callback =
    # New:       { bool user_feasible = true;
    #              if (mipsolver.submip) { expand via undoPrimal, check }
    #              else { check directly }
    #              if (!user_feasible) return false; }
    #            const bool execute_mip_solution_callback =
    string(REPLACE
      "  const bool execute_mip_solution_callback ="
      "  // Reject solutions that violate user lazy constraints (e.g. SEC/GSEC)\n  {\n    bool user_feasible;\n    if (mipsolver.submip) {\n      // Sub-MIP: sol is in reduced (presolved) space; expand to original\n      HighsSolution temp_sol;\n      temp_sol.col_value = sol;\n      temp_sol.value_valid = true;\n      postSolveStack.undoPrimal(*mipsolver.options_mip_, temp_sol);\n      user_feasible = HighsUserSeparator::isFeasible(temp_sol.col_value);\n    } else {\n      user_feasible = HighsUserSeparator::isFeasible(sol);\n    }\n    if (!user_feasible) return false;\n  }\n  const bool execute_mip_solution_callback ="
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

# ── 4. Patch HighsSearch.cpp: add HighsUserPropagator call after reduced-cost fixing ──
file(READ "${MIP_DIR}/HighsSearch.cpp" CONTENT3)

string(FIND "${CONTENT3}" "HighsUserPropagator" _found3)
if(_found3 EQUAL -1)
    # Add include for HighsUserPropagator.h
    string(REPLACE
      "#include \"mip/HighsSearch.h\""
      "#include \"mip/HighsSearch.h\"\n#include \"mip/HighsUserPropagator.h\""
      CONTENT3 "${CONTENT3}")

    # Insert propagator call after HighsRedcostFixing::propagateRedCost
    # and before localdom.propagate()
    string(REPLACE
      "HighsRedcostFixing::propagateRedCost(mipsolver, localdom, *lp);\n            localdom.propagate();"
      "HighsRedcostFixing::propagateRedCost(mipsolver, localdom, *lp);\n            HighsUserPropagator::propagate(localdom, mipsolver, *lp);\n            localdom.propagate();"
      CONTENT3 "${CONTENT3}")

    file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${CONTENT3}")
    message(STATUS "Applied HighsUserPropagator patch to HighsSearch.cpp")
else()
    message(STATUS "HighsSearch.cpp propagator patch already applied, skipping")
endif()

# ── 4a. Patch HighsLpRelaxation.h: add kBranching origin + methods ──
file(READ "${MIP_DIR}/HighsLpRelaxation.h" LPRELAX_H)

string(FIND "${LPRELAX_H}" "kBranching" _found_lprelax_h)
if(_found_lprelax_h EQUAL -1)
    # Add kBranching to the Origin enum (match the full enum context to avoid
    # accidentally replacing kCutPool in factory methods)
    string(REPLACE
      "      kCutPool,\n    };"
      "      kCutPool,\n      kBranching,\n    };"
      LPRELAX_H "${LPRELAX_H}")

    # Add static factory method branching() after model() factory
    string(REPLACE
      "static LpRow model(HighsInt index) { return LpRow{kModel, index, 0}; }"
      "static LpRow model(HighsInt index) { return LpRow{kModel, index, 0}; }\n    static LpRow branching() { return LpRow{kBranching, -1, 0}; }"
      LPRELAX_H "${LPRELAX_H}")

    # Add method declarations after removeObsoleteRows (unique, appears once)
    string(REPLACE
      "void removeObsoleteRows(bool notifyPool = true);"
      "void removeObsoleteRows(bool notifyPool = true);\n\n  HighsInt addBranchingRow(HighsInt nnz, const HighsInt* inds,\n                            const double* vals, double lo, double hi);\n\n  void removeLastBranchingRow();"
      LPRELAX_H "${LPRELAX_H}")

    file(WRITE "${MIP_DIR}/HighsLpRelaxation.h" "${LPRELAX_H}")
    message(STATUS "Applied kBranching patch to HighsLpRelaxation.h")
else()
    message(STATUS "HighsLpRelaxation.h kBranching patch already applied, skipping")
endif()

# ── 4b. Patch HighsLpRelaxation.cpp: implement addBranchingRow/removeLastBranchingRow ──
file(READ "${MIP_DIR}/HighsLpRelaxation.cpp" LPRELAX_CPP)

string(FIND "${LPRELAX_CPP}" "addBranchingRow" _found_lprelax_cpp)
if(_found_lprelax_cpp EQUAL -1)
    # Handle kBranching in get(): return len=0 (branching rows have no cut/model data)
    string(REPLACE
      "  switch (origin) {\n    case kCutPool:\n      mipsolver.mipdata_->cutpool.getCut(index, len, inds, vals);\n      break;\n    case kModel:\n      mipsolver.mipdata_->getRow(index, len, inds, vals);\n  };"
      "  switch (origin) {\n    case kCutPool:\n      mipsolver.mipdata_->cutpool.getCut(index, len, inds, vals);\n      break;\n    case kModel:\n      mipsolver.mipdata_->getRow(index, len, inds, vals);\n      break;\n    case kBranching:\n      len = 0; inds = nullptr; vals = nullptr;\n      break;\n  };"
      LPRELAX_CPP "${LPRELAX_CPP}")

    # Handle kBranching in getRowLen(): return 0
    string(REPLACE
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.getRowLength(index);\n    case kModel:\n      return mipsolver.mipdata_->ARstart_[index + 1] -\n             mipsolver.mipdata_->ARstart_[index];\n  };"
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.getRowLength(index);\n    case kModel:\n      return mipsolver.mipdata_->ARstart_[index + 1] -\n             mipsolver.mipdata_->ARstart_[index];\n    case kBranching:\n      return 0;\n  };"
      LPRELAX_CPP "${LPRELAX_CPP}")

    # Handle kBranching in isIntegral(): return false
    string(REPLACE
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.cutIsIntegral(index);\n    case kModel:\n      return (mipsolver.mipdata_->rowintegral[index] != 0);\n  };"
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.cutIsIntegral(index);\n    case kModel:\n      return (mipsolver.mipdata_->rowintegral[index] != 0);\n    case kBranching:\n      return false;\n  };"
      LPRELAX_CPP "${LPRELAX_CPP}")

    # Handle kBranching in getMaxAbsVal(): return 0
    string(REPLACE
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.getMaxAbsCutCoef(index);\n    case kModel:\n      return mipsolver.mipdata_->maxAbsRowCoef[index];\n  };"
      "  switch (origin) {\n    case kCutPool:\n      return mipsolver.mipdata_->cutpool.getMaxAbsCutCoef(index);\n    case kModel:\n      return mipsolver.mipdata_->maxAbsRowCoef[index];\n    case kBranching:\n      return 0.0;\n  };"
      LPRELAX_CPP "${LPRELAX_CPP}")

    # Handle kBranching in slackLower(): return row lower bound
    string(REPLACE
      "  switch (lprows[row].origin) {\n    case LpRow::kCutPool:\n      return mipsolver.mipdata_->domain.getMinCutActivity(\n          mipsolver.mipdata_->cutpool, lprows[row].index);\n    case LpRow::kModel:"
      "  switch (lprows[row].origin) {\n    case LpRow::kBranching:\n      return rowLower(row);\n    case LpRow::kCutPool:\n      return mipsolver.mipdata_->domain.getMinCutActivity(\n          mipsolver.mipdata_->cutpool, lprows[row].index);\n    case LpRow::kModel:"
      LPRELAX_CPP "${LPRELAX_CPP}")

    # Handle kBranching in slackUpper(): return rowupper
    string(REPLACE
      "  switch (lprows[row].origin) {\n    case LpRow::kCutPool:\n      return rowupper;\n    case LpRow::kModel:"
      "  switch (lprows[row].origin) {\n    case LpRow::kBranching:\n      return rowupper;\n    case LpRow::kCutPool:\n      return rowupper;\n    case LpRow::kModel:"
      LPRELAX_CPP "${LPRELAX_CPP}")

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
else()
    message(STATUS "HighsLpRelaxation.cpp branching row patch already applied, skipping")
endif()

# ── 5a. Patch HighsSearch.h: add hp_stack_pos and hp_data_idx to NodeData ──
file(READ "${MIP_DIR}/HighsSearch.h" SEARCH_H)

string(FIND "${SEARCH_H}" "hp_stack_pos" _found_search_h)
if(_found_search_h EQUAL -1)
    # Add fields to NodeData struct, after domgchgStackPos
    string(REPLACE
      "HighsInt domgchgStackPos;"
      "HighsInt domgchgStackPos;\n    HighsInt hp_stack_pos = 0;\n    HighsInt hp_data_idx = -1;"
      SEARCH_H "${SEARCH_H}")

    # Add include for HighsUserSeparator
    string(REPLACE
      "#include \"mip/HighsSearch.h\""
      "#include \"mip/HighsSearch.h\""
      SEARCH_H "${SEARCH_H}")

    file(WRITE "${MIP_DIR}/HighsSearch.h" "${SEARCH_H}")
    message(STATUS "Applied hp_stack_pos/hp_data_idx patch to HighsSearch.h")
else()
    message(STATUS "HighsSearch.h hyperplane fields already applied, skipping")
endif()

# ── 5b. Patch HighsSearch.cpp: hyperplane branching in branch() ──
file(READ "${MIP_DIR}/HighsSearch.cpp" SEARCH_CONTENT)

string(FIND "${SEARCH_CONTENT}" "HighsUserSeparator" _found_search)
if(_found_search EQUAL -1)
    # Add includes
    string(REPLACE
      "#include \"mip/HighsSearch.h\""
      "#include \"mip/HighsSearch.h\"\n#include \"mip/HighsUserSeparator.h\""
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # Insert hyperplane candidate evaluation + dynamic row branching
    # before normal child creation.
    string(REPLACE
      "  // finally open a new node with the branching decision added\n  // and remember that we have one open subtree left\n  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  currnode.opensubtrees = 1;\n  nodestack.emplace_back(\n      std::max(childLb, currnode.lower_bound), currnode.estimate,\n      currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n  nodestack.back().domgchgStackPos = domchgPos;\n\n  return NodeResult::kBranched;"
      "  // Hyperplane branching: evaluate candidates vs variable branching\n  const auto& branchCb = HighsUserSeparator::getBranchingCallback();\n  if (branchCb && !mipsolver.submip) {\n    auto candidates = branchCb(mipsolver);\n    if (!candidates.empty()) {\n      // Use pseudocost estimates for var_score (comparable with HP scoring).\n      // SB-based childLb/other_child_lb are only meaningful when minreliable > 0.\n      double var_score = 0.0;\n      if (currnode.branchingdecision.column >= 0) {\n        double vd = pseudocost.getPseudocostDown(\n            currnode.branchingdecision.column, currnode.branching_point);\n        double vu = pseudocost.getPseudocostUp(\n            currnode.branchingdecision.column, currnode.branching_point);\n        var_score = std::min(vd, vu) * std::max(vd, vu);\n      }\n\n      int winner = HighsUserSeparator::evaluateCandidates(\n          candidates, var_score, mipsolver, (int)getCurrentDepth());\n      if (winner >= 0) {\n        auto& cand = candidates[winner];\n        const auto& sol = lp->getSolution().col_value;\n        double lhs = 0.0;\n        for (int j = 0; j < (int)cand.indices.size(); ++j)\n          lhs += cand.values[j] * sol[cand.indices[j]];\n\n        // Store hyperplane data for flip\n        HighsInt data_idx = HighsUserSeparator::allocHyperplaneData();\n        auto& data = HighsUserSeparator::getHyperplaneData(data_idx);\n        data.indices = cand.indices;\n        data.values = cand.values;\n        data.branch_value = lhs;\n        currnode.hp_data_idx = data_idx;\n\n        // Add branching row (down branch: lhs <= floor)\n        HighsInt hpPos = HighsUserSeparator::stackSize();\n        HighsUserSeparator::pushBranching(\n            mipsolver.mipdata_->lp, data,\n            -kHighsInf, std::floor(lhs));\n\n        // Push no-op domain marker so backtrack() has a kBranching sentinel\n        HighsInt domchgPos = localdom.getDomainChangeStack().size();\n        localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n        currnode.opensubtrees = 1;\n        nodestack.emplace_back(\n            currnode.lower_bound, currnode.estimate,\n            currnode.nodeBasis, nullptr);\n        nodestack.back().domgchgStackPos = domchgPos;\n        nodestack.back().hp_stack_pos = hpPos;\n\n        return NodeResult::kBranched;\n      }\n    }\n  }\n\n  // finally open a new node with the branching decision added\n  // and remember that we have one open subtree left\n  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  currnode.opensubtrees = 1;\n  nodestack.emplace_back(\n      std::max(childLb, currnode.lower_bound), currnode.estimate,\n      currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n  nodestack.back().domgchgStackPos = domchgPos;\n  nodestack.back().hp_stack_pos = HighsUserSeparator::stackSize();\n\n  return NodeResult::kBranched;"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5c. Patch backtrack pop_back + localdom.backtrack(): restore branching rows ──
    # Pattern in HiGHS v1.13.1 has #ifndef NDEBUG around the backtrack() return value.
    # We match the exact pattern: nodestack.pop_back() followed by #ifndef block.
    # This pattern appears 3 times (backtrack, backtrackPlunge, backtrackUntilDepth).
    string(REPLACE
      "      nodestack.pop_back();\n#ifndef NDEBUG\n      HighsDomainChange branchchg =\n#endif\n          localdom.backtrack();"
      "      HighsInt hpPos_pop = nodestack.back().hp_stack_pos;\n      nodestack.pop_back();\n#ifndef NDEBUG\n      HighsDomainChange branchchg =\n#endif\n          localdom.backtrack();\n      HighsUserSeparator::restoreBranching(mipsolver.mipdata_->lp, hpPos_pop);"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # backtrackUntilDepth has a slightly different pattern (4-space indent)
    string(REPLACE
      "    nodestack.pop_back();\n\n#ifndef NDEBUG\n    HighsDomainChange branchchg =\n#endif\n        localdom.backtrack();"
      "    HighsInt hpPos_pop = nodestack.back().hp_stack_pos;\n    nodestack.pop_back();\n\n#ifndef NDEBUG\n    HighsDomainChange branchchg =\n#endif\n        localdom.backtrack();\n    HighsUserSeparator::restoreBranching(mipsolver.mipdata_->lp, hpPos_pop);"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5d. Patch flip sections in backtrack() and backtrackPlunge() ──
    # When hp_data_idx >= 0, we skip the branchingdecision flip, push a no-op
    # marker, and add the up-branch hyperplane row.
    # This pattern appears twice (backtrack + backtrackPlunge) with 4-space indent.
    string(REPLACE
      "    size_t numChangedCols = localdom.getChangedCols().size();\n    bool passStabilizerToChildNode =\n        orbitsValidInChildNode(currnode.branchingdecision);\n    localdom.changeBound(currnode.branchingdecision);"
      "    size_t numChangedCols = localdom.getChangedCols().size();\n    HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n    if (currnode.hp_data_idx >= 0) {\n      // Hyperplane flip: push no-op marker + up-branch row\n      auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n      localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n      HighsUserSeparator::pushBranching(\n          mipsolver.mipdata_->lp, flipData,\n          std::ceil(flipData.branch_value), kHighsInf);\n    } else {\n      localdom.changeBound(currnode.branchingdecision);\n    }\n    bool passStabilizerToChildNode =\n        (currnode.hp_data_idx < 0) && orbitsValidInChildNode(currnode.branchingdecision);"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # backtrackPlunge has HighsInt instead of size_t for numChangedCols
    string(REPLACE
      "    HighsInt numChangedCols = localdom.getChangedCols().size();\n    bool passStabilizerToChildNode =\n        orbitsValidInChildNode(currnode.branchingdecision);\n    localdom.changeBound(currnode.branchingdecision);"
      "    HighsInt numChangedCols = localdom.getChangedCols().size();\n    HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n    if (currnode.hp_data_idx >= 0) {\n      auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n      localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n      HighsUserSeparator::pushBranching(\n          mipsolver.mipdata_->lp, flipData,\n          std::ceil(flipData.branch_value), kHighsInf);\n    } else {\n      localdom.changeBound(currnode.branchingdecision);\n    }\n    bool passStabilizerToChildNode =\n        (currnode.hp_data_idx < 0) && orbitsValidInChildNode(currnode.branchingdecision);"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5e. Patch flip in backtrackUntilDepth() ──
    # This one has 2-space indent.
    string(REPLACE
      "  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  nodestack.emplace_back(\n      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n\n  lp->flushDomain(localdom);\n  nodestack.back().domgchgStackPos = domchgPos;"
      "  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n  HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n  if (currnode.hp_data_idx >= 0) {\n    auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n    localdom.changeBound({localdom.col_lower_[currnode.branchingdecision.column],\n          currnode.branchingdecision.column, HighsBoundType::kLower});\n    HighsUserSeparator::pushBranching(\n        mipsolver.mipdata_->lp, flipData,\n        std::ceil(flipData.branch_value), kHighsInf);\n  } else {\n    localdom.changeBound(currnode.branchingdecision);\n  }\n  bool passStabilizerToChildNode =\n      (currnode.hp_data_idx < 0) && orbitsValidInChildNode(currnode.branchingdecision);\n  nodestack.emplace_back(\n      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n\n  lp->flushDomain(localdom);\n  nodestack.back().domgchgStackPos = domchgPos;\n  nodestack.back().hp_stack_pos = hpPos_flip;"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5f. Set hp_stack_pos on flip child in backtrack() and backtrackPlunge() ──
    # After the flip's emplace_back, set hp_stack_pos. The pattern ends with:
    #   lp->flushDomain(localdom);
    #   nodestack.back().domgchgStackPos = domchgPos;
    #   break;
    # This appears in both backtrack() and backtrackPlunge().
    string(REPLACE
      "    lp->flushDomain(localdom);\n    nodestack.back().domgchgStackPos = domchgPos;\n    break;"
      "    lp->flushDomain(localdom);\n    nodestack.back().domgchgStackPos = domchgPos;\n    nodestack.back().hp_stack_pos = hpPos_flip;\n    break;"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5g. Restore branching rows when flip child is pruned or queued ──
    # When the flip child is pruned (infeasible/cutoff) or sent to the node queue,
    # localdom.backtrack() is called but branching rows must also be removed.
    # The prune pattern is: "localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight)"
    # This appears in backtrack() and backtrackPlunge().
    string(REPLACE
      "      localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());\n      continue;"
      "      HighsUserSeparator::restoreBranching(mipsolver.mipdata_->lp, hpPos_flip);\n      localdom.backtrack();\n      localdom.clearChangedCols(numChangedCols);\n      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());\n      continue;"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5h. Prevent hyperplane-branched nodes from going to the node queue ──
    # The node queue stores only domain changes, not LP branching rows.
    # Also add a guard AFTER the ancestor score loop since it can override nodeToQueue.
    string(REPLACE
      "    bool nodeToQueue = nodelb > mipsolver.mipdata_->optimality_limit;\n    // we check if switching"
      "    bool nodeToQueue = nodelb > mipsolver.mipdata_->optimality_limit;\n    if (currnode.hp_data_idx >= 0) nodeToQueue = false;  // hyperplane nodes can't be queued\n    // we check if switching"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # Add a second guard right before the nodeToQueue block
    string(REPLACE
      "    if (nodeToQueue) {\n      // if (!mipsolver.submip) printf(\"node goes to queue\\n\");"
      "    if (currnode.hp_data_idx >= 0) nodeToQueue = false;  // guard after ancestor score check\n    if (nodeToQueue) {\n      // if (!mipsolver.submip) printf(\"node goes to queue\\n\");"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5i. Suppress debug asserts for hyperplane branches ──
    # The asserts check branchchg matches branchingdecision. For hyperplane branches
    # the domain change is a no-op marker, not the branchingdecision.
    # Guard with hp_data_idx check.

    # Pattern in backtrack() and backtrackPlunge() (6-space indent, multi-line assert)
    string(REPLACE
      "      assert(\n          (branchchg.boundtype == HighsBoundType::kLower &&\n           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n          (branchchg.boundtype == HighsBoundType::kUpper &&\n           branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n      assert(branchchg.boundtype ==\n             nodestack.back().branchingdecision.boundtype);\n      assert(branchchg.column == nodestack.back().branchingdecision.column);"
      "      if (nodestack.back().hp_data_idx < 0) {\n      assert(\n          (branchchg.boundtype == HighsBoundType::kLower &&\n           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n          (branchchg.boundtype == HighsBoundType::kUpper &&\n           branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n      assert(branchchg.boundtype ==\n             nodestack.back().branchingdecision.boundtype);\n      assert(branchchg.column == nodestack.back().branchingdecision.column);\n      }"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # Pattern in backtrackUntilDepth() (4-space indent, single-line asserts)
    string(REPLACE
      "    assert(\n        (branchchg.boundtype == HighsBoundType::kLower &&\n         branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n        (branchchg.boundtype == HighsBoundType::kUpper &&\n         branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);\n    assert(branchchg.column == nodestack.back().branchingdecision.column);"
      "    if (nodestack.back().hp_data_idx < 0) {\n    assert(\n        (branchchg.boundtype == HighsBoundType::kLower &&\n         branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||\n        (branchchg.boundtype == HighsBoundType::kUpper &&\n         branchchg.boundval <= nodestack.back().branchingdecision.boundval));\n    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);\n    assert(branchchg.column == nodestack.back().branchingdecision.column);\n    }"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${SEARCH_CONTENT}")
    message(STATUS "Applied hyperplane branching patch to HighsSearch.cpp")
else()
    message(STATUS "HighsSearch.cpp hyperplane patch already applied, skipping")
endif()
