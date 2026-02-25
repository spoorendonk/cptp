#include "mip/HighsUserPropagator.h"

// Thread-local callback storage — each thread gets its own copy,
// enabling concurrent Highs instances on different threads.
namespace {
thread_local HighsUserPropagator::Callback callback_;
}  // namespace

void HighsUserPropagator::setCallback(Callback cb) {
  callback_ = std::move(cb);
}

void HighsUserPropagator::clearCallback() {
  callback_ = nullptr;
}

bool HighsUserPropagator::hasCallback() {
  return static_cast<bool>(callback_);
}

void HighsUserPropagator::propagate(HighsDomain& domain,
                                    const HighsMipSolver& mipsolver,
                                    const HighsLpRelaxation& lp) {
  if (callback_) {
    callback_(domain, mipsolver, lp);
  }
}
