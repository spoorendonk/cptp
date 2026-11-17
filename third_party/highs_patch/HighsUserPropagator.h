#ifndef MIP_HIGHS_USER_PROPAGATOR_H_
#define MIP_HIGHS_USER_PROPAGATOR_H_

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsUserHooks.h"

/// User-defined domain propagator injected into HiGHS's search loop.
/// The setter fills the calling thread's pending hooks, which a MIP solve
/// captures at run() (same pattern as HighsUserSeparator).
///
/// Called after reduced-cost fixing in each node evaluation.
/// The callback can fix variables via domain.changeBound().
class HighsUserPropagator {
public:
    using Callback = HighsUserHooks::Propagator;

    /// Set this thread's pending propagator callback. Call before Highs::run().
    static void setCallback(Callback cb);

    /// Clear this thread's pending propagator callback.
    static void clearCallback();

    /// Whether this thread has a pending propagator callback.
    static bool hasCallback();

    /// Run the MIP solve's propagator callback (called from HighsSearch.cpp).
    static void propagate(HighsDomain& domain, const HighsMipSolver& mipsolver,
                          const HighsLpRelaxation& lp);
};

#endif  // MIP_HIGHS_USER_PROPAGATOR_H_
