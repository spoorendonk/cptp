#ifndef MIP_HIGHS_USER_PROPAGATOR_H_
#define MIP_HIGHS_USER_PROPAGATOR_H_

#include <functional>

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"

/// User-defined domain propagator injected into HiGHS's search loop.
/// Uses thread-local callback storage (same pattern as HighsUserSeparator).
///
/// Called after reduced-cost fixing in each node evaluation.
/// The callback can fix variables via domain.changeBound().
class HighsUserPropagator {
 public:
  using Callback = std::function<void(HighsDomain& domain,
                                      const HighsMipSolver& mipsolver,
                                      const HighsLpRelaxation& lp)>;

  /// Set the thread-local propagator callback. Call before Highs::run().
  static void setCallback(Callback cb);

  /// Clear the thread-local callback.
  static void clearCallback();

  /// Run the propagator callback (called from HighsSearch.cpp).
  static void propagate(HighsDomain& domain,
                        const HighsMipSolver& mipsolver,
                        const HighsLpRelaxation& lp);
};

#endif  // MIP_HIGHS_USER_PROPAGATOR_H_
