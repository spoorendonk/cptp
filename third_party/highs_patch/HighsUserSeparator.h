#ifndef MIP_HIGHS_USER_SEPARATOR_H_
#define MIP_HIGHS_USER_SEPARATOR_H_

#include <functional>
#include <vector>

#include "mip/HighsSeparator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsMipSolver.h"

/// User-defined separator injected into HiGHS's separation loop.
/// Uses thread-local callback storage so multiple Highs instances can
/// run concurrently on different threads without interfering.
///
/// Thread-local variables are defined in HighsUserSeparator.cpp.
class HighsUserSeparator : public HighsSeparator {
 public:
  using Callback = std::function<void(const HighsLpRelaxation& lpRelaxation,
                                      HighsCutPool& cutpool,
                                      const HighsMipSolver& mipsolver)>;

  /// Solution feasibility check: returns true if the solution is feasible
  /// (no violated lazy constraints). Called from addIncumbent() before
  /// accepting any integer solution (including heuristic solutions).
  using FeasibilityCheck = std::function<bool(const std::vector<double>& sol)>;

  HighsUserSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, "User") {}

  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& /*lpAggregator*/,
                          HighsTransformedLp& /*transLp*/,
                          HighsCutPool& cutpool) override;

  /// Set the thread-local user separator callback. Call before Highs::run().
  static void setCallback(Callback cb);

  /// Set the thread-local solution feasibility check. Call before Highs::run().
  static void setFeasibilityCheck(FeasibilityCheck cb);

  /// Clear all thread-local callbacks.
  static void clearCallback();

  /// Check if a separator callback is registered (for lazy constraint loop).
  static bool hasCallback();

  /// Check if a solution is feasible w.r.t. user lazy constraints.
  /// Returns true if no check is registered or if the solution passes.
  static bool isFeasible(const std::vector<double>& sol);
};

#endif  // MIP_HIGHS_USER_SEPARATOR_H_
