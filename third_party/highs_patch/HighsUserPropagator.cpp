#include "mip/HighsUserPropagator.h"

// Thread-local callback storage — each thread gets its own copy,
// enabling concurrent Highs instances on different threads.
// The name is unique among the patch's sources: on Apple, HiGHS builds its
// library as a unity build, which can put this file and HighsUserSeparator.cpp
// (whose anonymous namespace has its own `callback_`) into one unit.
namespace {
thread_local HighsUserPropagator::Callback propagator_callback_;
}  // namespace

void HighsUserPropagator::setCallback(Callback cb) {
    propagator_callback_ = std::move(cb);
}

void HighsUserPropagator::clearCallback() {
    propagator_callback_ = nullptr;
}

bool HighsUserPropagator::hasCallback() {
    return static_cast<bool>(propagator_callback_);
}

void HighsUserPropagator::propagate(HighsDomain& domain, const HighsMipSolver& mipsolver,
                                    const HighsLpRelaxation& lp) {
    if (propagator_callback_) {
        propagator_callback_(domain, mipsolver, lp);
    }
}
