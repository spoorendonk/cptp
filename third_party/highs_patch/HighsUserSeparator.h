#ifndef MIP_HIGHS_USER_SEPARATOR_H_
#define MIP_HIGHS_USER_SEPARATOR_H_

#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "mip/HighsSeparator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
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

  /// Hyperplane candidate: sparse row defined by column indices + coefficients.
  struct HyperplaneCandidate {
      std::vector<HighsInt> indices;   // column indices
      std::vector<double> values;      // coefficients
      std::string name;                // key for pseudocost tracking
  };

  /// Hyperplane data stored in pool for backtrack/flip.
  struct HyperplaneData {
      std::vector<HighsInt> indices;
      std::vector<double> values;
      double branch_value;  // LHS value at branch point
  };

  /// Online pseudocost tracker for hyperplane candidates, keyed by name.
  struct HyperplanePseudocost {
      double cost_down = 0.0;
      double cost_up = 0.0;
      int samples_down = 0;
      int samples_up = 0;

      void addObservation(double frac_val, double obj_improvement, bool is_up) {
          double delta = is_up ? (std::ceil(frac_val) - frac_val)
                               : (frac_val - std::floor(frac_val));
          if (delta < 1e-12) return;
          double unit_gain = obj_improvement / delta;
          if (is_up) {
              cost_up += (unit_gain - cost_up) / (++samples_up);
          } else {
              cost_down += (unit_gain - cost_down) / (++samples_down);
          }
      }

      double score(double frac_val) const {
          double frac = frac_val - std::floor(frac_val);
          double d = cost_down * frac;
          double u = cost_up * (1.0 - frac);
          return std::min(d, u) * std::max(d, u);
      }

      bool reliable(int min_reliable) const {
          return samples_up >= min_reliable && samples_down >= min_reliable;
      }
  };

  /// Callback: return candidate hyperplanes for branching evaluation.
  using BranchingCallback = std::function<std::vector<HyperplaneCandidate>(
      const HighsMipSolver& mipsolver)>;

  /// Configuration for strong branching evaluation.
  struct StrongBranchConfig {
      int max_depth = 0;           // SB depth limit (0 = no SB, dual-score only)
      int iter_limit = 100;        // simplex iterations per trial solve
      int min_reliable = 4;        // pseudocost samples needed before trusted
      int max_sb_candidates = 3;   // top-k from dual pre-score for SB
  };

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

  /// Set the global branching callback. Call before Highs::run().
  static void setBranchingCallback(BranchingCallback cb);

  /// Get the global branching callback.
  static const BranchingCallback& getBranchingCallback();

  /// Set strong branching configuration.
  static void setStrongBranchConfig(StrongBranchConfig cfg);

  /// Get strong branching configuration.
  static const StrongBranchConfig& getStrongBranchConfig();

  /// Check if a separator callback is registered (for lazy constraint loop).
  static bool hasCallback();

  /// Check if a solution is feasible w.r.t. user lazy constraints.
  /// Returns true if no check is registered or if the solution passes.
  static bool isFeasible(const std::vector<double>& sol);

  /// 3-tier evaluation pipeline: dual pre-score -> pseudocost -> sequential SB.
  /// Returns index of winning candidate, or -1 if variable branching wins.
  static int evaluateCandidates(
      std::vector<HyperplaneCandidate>& candidates,
      double var_score,
      const HighsMipSolver& mipsolver,
      int depth);

  // -- Hyperplane branching stack/pool management --

  /// Allocate a new HyperplaneData entry in the pool. Returns its index.
  static HighsInt allocHyperplaneData();

  /// Get a reference to pool entry by index.
  static HyperplaneData& getHyperplaneData(HighsInt idx);

  /// Current pool size (= next alloc index).
  static HighsInt poolSize();

  /// Current stack size (number of active branching rows).
  static HighsInt stackSize();

  /// Push a branching row: add row to LP + push stack marker.
  static void pushBranching(HighsLpRelaxation& lp,
                            const HyperplaneData& data,
                            double lo, double hi);

  /// Restore branching rows: remove rows until stack matches target_pos.
  static void restoreBranching(HighsLpRelaxation& lp, HighsInt target_pos);

  /// Clear stack and pool (call between solves).
  static void clearStack();

 private:
  // callback_ and feasibility_check_ are thread_local in HighsUserSeparator.cpp
  static BranchingCallback branching_callback_;
  static StrongBranchConfig strong_branch_config_;
  static std::mutex mutex_;
  static std::unordered_map<std::string, HyperplanePseudocost> hyperplane_pseudocosts_;

  // Hyperplane branching: row lifecycle tracking
  static std::vector<HighsInt> hp_stack_;        // tracks active branching rows
  static std::vector<HyperplaneData> hp_data_pool_;  // coefficients for flip
};

#endif  // MIP_HIGHS_USER_SEPARATOR_H_
