#include "mip/HighsUserPropagator.h"

#include <utility>

void HighsUserPropagator::setCallback(Callback cb) {
    HighsUserHooks::pending().propagator = std::move(cb);
}

void HighsUserPropagator::clearCallback() {
    HighsUserHooks::pending().propagator = nullptr;
}

bool HighsUserPropagator::hasCallback() {
    return static_cast<bool>(HighsUserHooks::pending().propagator);
}

void HighsUserPropagator::propagate(HighsDomain& domain, const HighsMipSolver& mipsolver,
                                    const HighsLpRelaxation& lp) {
    if (mipsolver.user_hooks_ && mipsolver.user_hooks_->propagator) {
        mipsolver.user_hooks_->propagator(domain, mipsolver, lp);
    }
}
