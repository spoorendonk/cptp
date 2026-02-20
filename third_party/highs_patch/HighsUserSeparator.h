#ifndef MIP_HIGHS_USER_SEPARATOR_H_
#define MIP_HIGHS_USER_SEPARATOR_H_

#include <functional>
#include <mutex>
#include <vector>

#include "mip/HighsSeparator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsMipSolver.h"

/// User-defined separator injected into HiGHS's separation loop.
/// Uses a singleton callback registry so the callback can be set from
/// outside HiGHS before calling run().
///
/// The static members are defined in HighsUserSeparator.cpp (part of libhighs)
/// to ensure a single definition across the entire program.
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

  /// Set the global user separator callback. Call before Highs::run().
  static void setCallback(Callback cb);

  /// Set the global solution feasibility check. Call before Highs::run().
  static void setFeasibilityCheck(FeasibilityCheck cb);

  /// Clear all global callbacks.
  static void clearCallback();

  /// Check if a separator callback is registered (for lazy constraint loop).
  static bool hasCallback();

  /// Check if a solution is feasible w.r.t. user lazy constraints.
  /// Returns true if no check is registered or if the solution passes.
  static bool isFeasible(const std::vector<double>& sol);

 private:
  static Callback callback_;
  static FeasibilityCheck feasibility_check_;
  static std::mutex mutex_;
};

#endif  // MIP_HIGHS_USER_SEPARATOR_H_
