#include "mip/HighsUserSeparator.h"

// Thread-local callback storage — each thread gets its own copy,
// enabling concurrent Highs instances on different threads.
namespace {
thread_local HighsUserSeparator::Callback callback_;
thread_local HighsUserSeparator::FeasibilityCheck feasibility_check_;
}  // namespace

void HighsUserSeparator::separateLpSolution(
    HighsLpRelaxation& lpRelaxation,
    HighsLpAggregator& /*lpAggregator*/,
    HighsTransformedLp& /*transLp*/,
    HighsCutPool& cutpool) {
  if (callback_) {
    callback_(lpRelaxation, cutpool, lpRelaxation.getMipSolver());
  }
}

void HighsUserSeparator::setCallback(Callback cb) {
  callback_ = std::move(cb);
}

void HighsUserSeparator::setFeasibilityCheck(FeasibilityCheck cb) {
  feasibility_check_ = std::move(cb);
}

void HighsUserSeparator::clearCallback() {
  callback_ = nullptr;
  feasibility_check_ = nullptr;
}

bool HighsUserSeparator::hasCallback() {
  return static_cast<bool>(callback_);
}

bool HighsUserSeparator::isFeasible(const std::vector<double>& sol) {
  if (!feasibility_check_) return true;
  return feasibility_check_(sol);
}
