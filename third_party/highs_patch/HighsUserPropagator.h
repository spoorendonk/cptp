#ifndef MIP_HIGHS_USER_PROPAGATOR_H_
#define MIP_HIGHS_USER_PROPAGATOR_H_

#include <functional>
#include <mutex>

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"

/// User-defined domain propagator injected into HiGHS's search loop.
/// Uses a singleton callback registry (same pattern as HighsUserSeparator).
///
/// Called after reduced-cost fixing in each node evaluation.
/// The callback can fix variables via domain.changeBound().
class HighsUserPropagator {
 public:
  using Callback = std::function<void(HighsDomain& domain,
                                      const HighsMipSolver& mipsolver,
                                      const HighsLpRelaxation& lp)>;

  /// Set the global propagator callback. Call before Highs::run().
  static void setCallback(Callback cb);

  /// Clear the global callback.
  static void clearCallback();

  /// Run the propagator callback (called from HighsSearch.cpp).
  static void propagate(HighsDomain& domain,
                        const HighsMipSolver& mipsolver,
                        const HighsLpRelaxation& lp);

 private:
  static Callback callback_;
  static std::mutex mutex_;
};

#endif  // MIP_HIGHS_USER_PROPAGATOR_H_
