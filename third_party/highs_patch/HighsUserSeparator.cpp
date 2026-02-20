#include "mip/HighsUserSeparator.h"

// Static member definitions — single copy linked into libhighs.
HighsUserSeparator::Callback HighsUserSeparator::callback_;
HighsUserSeparator::FeasibilityCheck HighsUserSeparator::feasibility_check_;
std::mutex HighsUserSeparator::mutex_;

void HighsUserSeparator::separateLpSolution(
    HighsLpRelaxation& lpRelaxation,
    HighsLpAggregator& /*lpAggregator*/,
    HighsTransformedLp& /*transLp*/,
    HighsCutPool& cutpool) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (callback_) {
    callback_(lpRelaxation, cutpool, lpRelaxation.getMipSolver());
  }
}

void HighsUserSeparator::setCallback(Callback cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = std::move(cb);
}

void HighsUserSeparator::setFeasibilityCheck(FeasibilityCheck cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  feasibility_check_ = std::move(cb);
}

void HighsUserSeparator::clearCallback() {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = nullptr;
  feasibility_check_ = nullptr;
}

bool HighsUserSeparator::hasCallback() {
  std::lock_guard<std::mutex> lock(mutex_);
  return static_cast<bool>(callback_);
}

bool HighsUserSeparator::isFeasible(const std::vector<double>& sol) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!feasibility_check_) return true;
  return feasibility_check_(sol);
}
