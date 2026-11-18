#ifndef MIP_HIGHS_USER_SEPARATOR_H_
#define MIP_HIGHS_USER_SEPARATOR_H_

#include "mip/HighsCutPool.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsSeparator.h"
#include "mip/HighsUserHooks.h"

#include <vector>

/// User-defined separator injected into HiGHS's separation loop.
///
/// The setters fill the calling thread's pending hooks (HighsUserHooks), so
/// several Highs instances can run concurrently on different threads without
/// interfering. A MIP solve reads the hooks it captured at run()
/// (HighsMipSolver::user_hooks_), on whatever thread its tree workers run.
class HighsUserSeparator : public HighsSeparator {
public:
    using Callback = HighsUserHooks::Separator;

    /// Solution feasibility check: returns true if the solution is feasible
    /// (no violated lazy constraints). Called from addIncumbent() before
    /// accepting any integer solution (including heuristic solutions), on the
    /// main solver and on every worker of the parallel tree search.
    using FeasibilityCheck = HighsUserHooks::FeasibilityCheck;

    HighsUserSeparator(const HighsMipSolver& mipsolver) : HighsSeparator(mipsolver, "User") {}

    void separateLpSolution(HighsLpRelaxation& lpRelaxation, HighsLpAggregator& /*lpAggregator*/,
                            HighsTransformedLp& /*transLp*/, HighsCutPool& cutpool) override;

    /// Set this thread's pending separator callback. Call before Highs::run().
    static void setCallback(Callback cb);

    /// Set this thread's pending solution feasibility check. Call before Highs::run().
    static void setFeasibilityCheck(FeasibilityCheck cb);

    /// Set this thread's pending root cut-loop limits: the round cap and the
    /// stall limit, each HighsUserHooks::kRootLimitHighs (HiGHS's own),
    /// kRootLimitUnlimited or a positive value. Call before Highs::run().
    static void setRootSeparationLimits(int max_rounds, int stall_limit);

    /// Clear this thread's pending separator callback and feasibility check,
    /// and reset its root cut-loop limits to HiGHS's.
    static void clearCallback();

    /// The round cap of the MIP solve's root cut loop, given HiGHS's own:
    /// the solve's root limit (kHighsIInf when unlimited), or highs_rounds
    /// for a sub-MIP or when the solve keeps HiGHS's.
    static HighsInt rootSepaRounds(const HighsMipSolver& mipsolver, HighsInt highs_rounds);

    /// The stall limit of the MIP solve's root cut loop, as rootSepaRounds.
    static HighsInt rootSepaStall(const HighsMipSolver& mipsolver, HighsInt highs_stall);

    /// Whether this thread has a pending separator callback.
    static bool hasCallback();

    /// Check a solution against this thread's pending feasibility check.
    /// Returns true if no check is pending or if the solution passes.
    static bool isFeasible(const std::vector<double>& sol);

    /// Whether the MIP solve has a separator callback (for the lazy
    /// constraint loop on integral LP solutions).
    static bool hasCallback(const HighsMipSolver& mipsolver);

    /// Check a solution of the MIP solve against its lazy constraints. A
    /// sub-MIP's solution is in its presolved space and is expanded to the
    /// original space first. Returns true if the solve has no feasibility
    /// check or if the solution passes.
    static bool isFeasible(const HighsMipSolver& mipsolver, const std::vector<double>& sol);

    /// Expand a solution of the MIP solve to the top-level model's columns,
    /// undoing the presolve of every sub-MIP level (HiGHS nests sub-MIPs).
    static std::vector<double> toOriginalSpace(const HighsMipSolver& mipsolver,
                                               const std::vector<double>& sol);

    /// The top-level model's column of each of the solve's first num_col
    /// columns, through every sub-MIP level.
    static std::vector<HighsInt> originalColumns(const HighsMipSolver& mipsolver, HighsInt num_col);

    /// Record that addIncumbent() rejected a solution because it violates the
    /// user lazy constraints. thread_local: the mark and its take happen in one
    /// synchronous call on one thread (HighsSearch::evaluateNode), on the main
    /// solver or on a tree worker.
    static void markIncumbentRejected();

    /// Return whether a rejection was recorded since the last call, and clear it.
    static bool takeIncumbentRejected();
};

#endif  // MIP_HIGHS_USER_SEPARATOR_H_
