#include "mip/HighsUserPropagator.h"

// Static member definitions — single copy linked into libhighs.
HighsUserPropagator::Callback HighsUserPropagator::callback_;
std::mutex HighsUserPropagator::mutex_;

void HighsUserPropagator::setCallback(Callback cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = std::move(cb);
}

void HighsUserPropagator::clearCallback() {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = nullptr;
}

void HighsUserPropagator::propagate(HighsDomain& domain,
                                    const HighsMipSolver& mipsolver,
                                    const HighsLpRelaxation& lp) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (callback_) {
    callback_(domain, mipsolver, lp);
  }
}
