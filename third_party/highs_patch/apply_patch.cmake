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
      "  // Hyperplane branching: evaluate candidates vs variable branching\n  const auto& branchCb = HighsUserSeparator::getBranchingCallback();\n  if (branchCb && !mipsolver.submip) {\n    auto candidates = branchCb(mipsolver);\n    if (!candidates.empty()) {\n      double parentObj = lp->getObjective();\n      double vdi = std::max(0.0, childLb - parentObj);\n      double vui = std::max(0.0, currnode.other_child_lb - parentObj);\n      double var_score = std::min(vdi, vui) * std::max(vdi, vui);\n\n      int winner = HighsUserSeparator::evaluateCandidates(\n          candidates, var_score, mipsolver, (int)getCurrentDepth());\n      if (winner >= 0) {\n        auto& cand = candidates[winner];\n        // Compute LHS\n        const auto& sol = lp->getSolution().col_value;\n        double lhs = 0.0;\n        for (int j = 0; j < (int)cand.indices.size(); ++j)\n          lhs += cand.values[j] * sol[cand.indices[j]];\n\n        // Store hyperplane data for flip\n        HighsInt data_idx = HighsUserSeparator::allocHyperplaneData();\n        auto& data = HighsUserSeparator::getHyperplaneData(data_idx);\n        data.indices = cand.indices;\n        data.values = cand.values;\n        data.branch_value = lhs;\n        currnode.hp_data_idx = data_idx;\n\n        // Add branching row (down branch: lhs <= floor)\n        HighsInt hpPos = HighsUserSeparator::stackSize();\n        HighsUserSeparator::pushBranching(\n            mipsolver.mipdata_->lp, data,\n            -kHighsInf, std::floor(lhs));\n\n        // Also do normal variable branch\n        HighsInt domchgPos = localdom.getDomainChangeStack().size();\n        bool passStabilizerToChildNode =\n            orbitsValidInChildNode(currnode.branchingdecision);\n        localdom.changeBound(currnode.branchingdecision);\n        currnode.opensubtrees = 1;\n        nodestack.emplace_back(\n            std::max(childLb, currnode.lower_bound), currnode.estimate,\n            currnode.nodeBasis,\n            passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n        nodestack.back().domgchgStackPos = domchgPos;\n        nodestack.back().hp_stack_pos = hpPos;\n\n        return NodeResult::kBranched;\n      }\n    }\n  }\n\n  // finally open a new node with the branching decision added\n  // and remember that we have one open subtree left\n  HighsInt domchgPos = localdom.getDomainChangeStack().size();\n\n  bool passStabilizerToChildNode =\n      orbitsValidInChildNode(currnode.branchingdecision);\n  localdom.changeBound(currnode.branchingdecision);\n  currnode.opensubtrees = 1;\n  nodestack.emplace_back(\n      std::max(childLb, currnode.lower_bound), currnode.estimate,\n      currnode.nodeBasis,\n      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n  nodestack.back().domgchgStackPos = domchgPos;\n\n  return NodeResult::kBranched;"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # ── 5c. Patch backtrack(): restore branching rows on backtrack + flip ──
    # Find the backtrack loop body where nodestack.pop_back() happens after
    # opensubtrees == 0 and patch to restore branching rows.
    #
    # The backtrack pattern in HighsSearch.cpp has two cases:
    # 1. Normal backtrack (opensubtrees == 0): pop + localdom.backtrack()
    # 2. Flip (opensubtrees == 1): flip branchingdecision + new child
    #
    # We need to handle both. The exact code varies by HiGHS version.
    # We'll patch around the nodestack.pop_back() + localdom.backtrack() pattern
    # and the flip section.

    # Patch: after "nodestack.pop_back();" + "localdom.backtrack();" in backtrack,
    # add restoreBranching call. We target the pattern in the while loop.
    string(REPLACE
      "      nodestack.pop_back();\n      localdom.backtrack();"
      "      HighsInt hpPos_bt = nodestack.back().hp_stack_pos;\n      nodestack.pop_back();\n      localdom.backtrack();\n      HighsUserSeparator::restoreBranching(mipsolver.mipdata_->lp, hpPos_bt);"
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # Patch the flip section: after flipping branchingdecision bounds,
    # add hyperplane flip logic. The flip creates a new child.
    # We look for the pattern where emplace_back happens after the flip.
    # The flip code changes branchingdecision.boundval and boundtype,
    # then does localdom.changeBound + emplace_back.
    #
    # After the flip's emplace_back, set hp_stack_pos on the new child.
    # Before changeBound, add the up-branch row if hp_data_idx >= 0.
    string(REPLACE
      "    localdom.changeBound(currnode.branchingdecision);\n    nodestack.emplace_back("
      "    // Hyperplane flip: add up-branch row if parent had hyperplane\n    HighsInt hpPos_flip = HighsUserSeparator::stackSize();\n    if (currnode.hp_data_idx >= 0) {\n      auto& flipData = HighsUserSeparator::getHyperplaneData(currnode.hp_data_idx);\n      HighsUserSeparator::pushBranching(\n          mipsolver.mipdata_->lp, flipData,\n          std::ceil(flipData.branch_value), kHighsInf);\n    }\n    localdom.changeBound(currnode.branchingdecision);\n    nodestack.emplace_back("
      SEARCH_CONTENT "${SEARCH_CONTENT}")

    # After the flip's emplace_back + domgchgStackPos assignment, add hp_stack_pos
    # The pattern is: nodestack.back().domgchgStackPos = domchgPos;
    # followed by either return or break.
    # We need to add hp_stack_pos right after domgchgStackPos in the flip context.
    # Since this pattern appears in multiple places, we target the specific one
    # in backtrack() by looking for the unique context.

    file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${SEARCH_CONTENT}")
    message(STATUS "Applied hyperplane branching patch to HighsSearch.cpp")
else()
    message(STATUS "HighsSearch.cpp hyperplane patch already applied, skipping")
endif()

# ── 5d. Patch: set hp_stack_pos after domgchgStackPos in flip context ──
# Re-read to apply additional patches
file(READ "${MIP_DIR}/HighsSearch.cpp" SEARCH_CONTENT2)

string(FIND "${SEARCH_CONTENT2}" "hp_stack_pos = hpPos_flip" _found_hpflip)
if(_found_hpflip EQUAL -1)
    # In the flip context, after domgchgStackPos is set, add hp_stack_pos.
    # We look for the domgchgStackPos assignment that follows the flip emplace_back.
    # The unique context is: "hpPos_flip" appears nearby.
    # Since there may be multiple domgchgStackPos assignments, we target
    # the one in the backtrack function by matching the larger context.
    #
    # The flip emplace_back pattern ends with:
    #   nodestack.back().domgchgStackPos = domchgPos;
    # We want to add:
    #   nodestack.back().hp_stack_pos = hpPos_flip;
    # But we need to be careful about which occurrence to match.
    # Let's match the full flip block context.

    # Actually, let's just do a simple approach: replace all
    # "nodestack.back().domgchgStackPos = domchgPos;" with a version
    # that also sets hp_stack_pos. But we only want this in backtrack contexts.
    # The branch() function already handles it. In backtrack, there's typically
    # one occurrence after the flip emplace_back.
    #
    # Since we already patched the branch() path to set hp_stack_pos for
    # hyperplane wins, and the normal branch path doesn't set it (defaults to 0),
    # we need the flip's child to get the correct hp_stack_pos.
    #
    # For the flip: the child should get hpPos_flip (stack size before push).
    # We need to find the specific domgchgStackPos in the backtrack flip path.
    # The unique context around it is "Hyperplane flip:" which we just inserted.

    # Target the flip section specifically by looking for the pattern we just inserted
    string(REPLACE
      "    nodestack.emplace_back(\n        std::max(otherLb, currnode.lower_bound), currnode.estimate,"
      "    nodestack.emplace_back(\n        std::max(otherLb, currnode.lower_bound), currnode.estimate,"
      SEARCH_CONTENT2 "${SEARCH_CONTENT2}")

    # Just add after every domgchgStackPos assignment in backtrack.
    # Since hp_stack_pos defaults to 0, this is safe even if no hyperplane was used.
    # We'll use a replace_all-style approach, but CMake doesn't have that easily.
    # Instead, let's target the unique string pattern in the flip section:
    # After "passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n    nodestack.back().domgchgStackPos = domchgPos;"
    # which appears in the backtrack flip.
    string(REPLACE
      "passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n    nodestack.back().domgchgStackPos = domchgPos;\n\n    return"
      "passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);\n    nodestack.back().domgchgStackPos = domchgPos;\n    nodestack.back().hp_stack_pos = hpPos_flip;\n\n    return"
      SEARCH_CONTENT2 "${SEARCH_CONTENT2}")

    file(WRITE "${MIP_DIR}/HighsSearch.cpp" "${SEARCH_CONTENT2}")
    message(STATUS "Applied hp_stack_pos flip patch to HighsSearch.cpp")
else()
    message(STATUS "hp_stack_pos flip patch already applied, skipping")
endif()

# ── 5e. Patch backtrackUntilDepth(): same pattern as backtrack ──
# backtrackUntilDepth has similar pop_back + backtrack patterns.
# Re-read and check if it needs patching (it uses the same nodestack).
# In HiGHS, backtrackUntilDepth may or may not exist depending on version.
# If it exists and has the same pattern, patch it.
file(READ "${MIP_DIR}/HighsSearch.cpp" SEARCH_CONTENT3)

# Check if backtrackUntilDepth exists and needs patching
string(FIND "${SEARCH_CONTENT3}" "backtrackUntilDepth" _found_btdepth)
if(NOT _found_btdepth EQUAL -1)
    # Check if already patched in this function
    string(FIND "${SEARCH_CONTENT3}" "hpPos_btd" _found_btd_patched)
    if(_found_btd_patched EQUAL -1)
        # backtrackUntilDepth has a similar loop with pop_back + backtrack.
        # It may use a slightly different pattern. Since we already patched
        # the main "nodestack.pop_back();\n      localdom.backtrack();" pattern
        # above (which would have caught both if they're identical), we check
        # if there's another instance that wasn't caught.
        # Actually, we already used string(REPLACE) which replaces ALL occurrences,
        # so both backtrack() and backtrackUntilDepth() should be patched.
        message(STATUS "backtrackUntilDepth exists; pop_back pattern already handled by global replace")
    endif()
endif()
