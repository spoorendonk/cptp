#include "mip/HighsUserSeparator.h"

#include "mip/HighsMipSolverData.h"

#include <numeric>
#include <utility>

namespace {
// The calling thread's pending hooks, captured by the next top-level
// HighsMipSolver::run() on this thread.
thread_local HighsUserHooks pending_hooks_;
thread_local bool incumbent_rejected_ = false;
}  // namespace

HighsUserHooks& HighsUserHooks::pending() {
    return pending_hooks_;
}

std::shared_ptr<const HighsUserHooks> HighsUserHooks::capture() {
    const HighsUserHooks& hooks = pending_hooks_;
    if (!hooks.separator && !hooks.feasibility_check && !hooks.propagator &&
        hooks.root_sepa_rounds == HighsUserHooks::kRootLimitHighs &&
        hooks.root_sepa_stall == HighsUserHooks::kRootLimitHighs)
        return nullptr;
    return std::make_shared<const HighsUserHooks>(hooks);
}

void HighsUserSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& /*lpAggregator*/,
                                            HighsTransformedLp& /*transLp*/,
                                            HighsCutPool& cutpool) {
    const HighsMipSolver& mipsolver = lpRelaxation.getMipSolver();
    if (hasCallback(mipsolver)) {
        mipsolver.user_hooks_->separator(lpRelaxation, cutpool, mipsolver);
    }
}

void HighsUserSeparator::setCallback(Callback cb) {
    pending_hooks_.separator = std::move(cb);
}

void HighsUserSeparator::setFeasibilityCheck(FeasibilityCheck cb) {
    pending_hooks_.feasibility_check = std::move(cb);
}

void HighsUserSeparator::setRootSeparationLimits(int max_rounds, int stall_limit) {
    pending_hooks_.root_sepa_rounds = max_rounds;
    pending_hooks_.root_sepa_stall = stall_limit;
}

void HighsUserSeparator::clearCallback() {
    pending_hooks_.separator = nullptr;
    pending_hooks_.feasibility_check = nullptr;
    pending_hooks_.root_sepa_rounds = HighsUserHooks::kRootLimitHighs;
    pending_hooks_.root_sepa_stall = HighsUserHooks::kRootLimitHighs;
}

namespace {
// A root limit of the solve's hooks in place of HiGHS's own (HighsUserHooks).
HighsInt rootLimit(const HighsMipSolver& mipsolver, int HighsUserHooks::* limit,
                   HighsInt highs_value) {
    if (mipsolver.submip || !mipsolver.user_hooks_)
        return highs_value;
    const int value = (*mipsolver.user_hooks_).*limit;
    if (value == HighsUserHooks::kRootLimitUnlimited)
        return kHighsIInf;
    if (value > 0)
        return value;
    return highs_value;
}
}  // namespace

HighsInt HighsUserSeparator::rootSepaRounds(const HighsMipSolver& mipsolver,
                                            HighsInt highs_rounds) {
    return rootLimit(mipsolver, &HighsUserHooks::root_sepa_rounds, highs_rounds);
}

HighsInt HighsUserSeparator::rootSepaStall(const HighsMipSolver& mipsolver, HighsInt highs_stall) {
    return rootLimit(mipsolver, &HighsUserHooks::root_sepa_stall, highs_stall);
}

bool HighsUserSeparator::hasCallback() {
    return static_cast<bool>(pending_hooks_.separator);
}

bool HighsUserSeparator::isFeasible(const std::vector<double>& sol) {
    if (!pending_hooks_.feasibility_check)
        return true;
    return pending_hooks_.feasibility_check(sol);
}

bool HighsUserSeparator::hasCallback(const HighsMipSolver& mipsolver) {
    return mipsolver.user_hooks_ && mipsolver.user_hooks_->separator;
}

bool HighsUserSeparator::isFeasible(const HighsMipSolver& mipsolver,
                                    const std::vector<double>& sol) {
    if (!mipsolver.user_hooks_ || !mipsolver.user_hooks_->feasibility_check)
        return true;
    const auto& check = mipsolver.user_hooks_->feasibility_check;
    if (!mipsolver.submip)
        return check(sol);
    return check(toOriginalSpace(mipsolver, sol));
}

std::vector<double> HighsUserSeparator::toOriginalSpace(const HighsMipSolver& mipsolver,
                                                        const std::vector<double>& sol) {
    std::vector<double> col_value = sol;
    // Each sub-MIP's model is its parent's presolved model, so one undo per
    // level. The walk stops below the top-level solve, whose presolve is off;
    // a sub-MIP runs on one thread, inside its parent's call, so the stacks
    // read here are not in use elsewhere.
    for (const HighsMipSolver* level = &mipsolver; level->submip; level = level->parent_mip_) {
        // A fresh solution per level: undo sizes row_value for one level only.
        HighsSolution expanded;
        expanded.col_value = std::move(col_value);
        expanded.value_valid = true;
        level->mipdata_->postSolveStack.undoPrimal(*level->options_mip_, expanded);
        col_value = std::move(expanded.col_value);
    }
    return col_value;
}

std::vector<HighsInt> HighsUserSeparator::originalColumns(const HighsMipSolver& mipsolver,
                                                          HighsInt num_col) {
    std::vector<HighsInt> cols(static_cast<size_t>(num_col));
    std::iota(cols.begin(), cols.end(), HighsInt{0});
    for (const HighsMipSolver* level = &mipsolver; level->submip; level = level->parent_mip_) {
        const auto& stack = level->mipdata_->postSolveStack;
        for (HighsInt& col : cols)
            col = stack.getOrigColIndex(col);
    }
    return cols;
}

void HighsUserSeparator::markIncumbentRejected() {
    incumbent_rejected_ = true;
}

bool HighsUserSeparator::takeIncumbentRejected() {
    const bool rejected = incumbent_rejected_;
    incumbent_rejected_ = false;
    return rejected;
}
