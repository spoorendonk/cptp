#include "mip/HighsUserSeparator.h"

#include "mip/HighsMipSolverData.h"

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
    if (!hooks.separator && !hooks.feasibility_check && !hooks.propagator)
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

void HighsUserSeparator::clearCallback() {
    pending_hooks_.separator = nullptr;
    pending_hooks_.feasibility_check = nullptr;
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
    // Sub-MIP: sol is in the reduced (presolved) space; expand to original.
    HighsSolution temp_sol;
    temp_sol.col_value = sol;
    temp_sol.value_valid = true;
    mipsolver.mipdata_->postSolveStack.undoPrimal(*mipsolver.options_mip_, temp_sol);
    return check(temp_sol.col_value);
}

void HighsUserSeparator::markIncumbentRejected() {
    incumbent_rejected_ = true;
}

bool HighsUserSeparator::takeIncumbentRejected() {
    const bool rejected = incumbent_rejected_;
    incumbent_rejected_ = false;
    return rejected;
}
