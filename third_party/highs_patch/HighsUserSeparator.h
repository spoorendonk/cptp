#ifndef MIP_HIGHS_USER_SEPARATOR_H_
#define MIP_HIGHS_USER_SEPARATOR_H_

#include "mip/HighsCutPool.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsSeparator.h"
#include "mip/HighsUserHooks.h"

#include <thread>
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

    /// Set this thread's pending root task observer (HighsUserHooks). Call
    /// before Highs::run().
    static void setRootTaskObserver(HighsUserHooks::RootTaskObserver observer);

    /// Clear this thread's pending separator callback, feasibility check and
    /// root task observer, and reset its root cut-loop limits to HiGHS's.
    static void clearCallback();

    /// Reports one root task (HighsUserHooks::root_task_observer) for its
    /// lifetime, when it runs on a thread other than `owner`, the one that
    /// spawned it: a task the owner runs itself at the sync keeps no other
    /// thread busy. Lives in the task's body, since a task cancelled before
    /// it starts never runs it.
    class RootTaskScope {
    public:
        RootTaskScope(const HighsMipSolver& mipsolver, std::thread::id owner);
        ~RootTaskScope();
        RootTaskScope(const RootTaskScope&) = delete;
        RootTaskScope& operator=(const RootTaskScope&) = delete;

    private:
        const HighsUserHooks::RootTaskObserver* observer_ = nullptr;
    };

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

    /// How the expanded value of a top-level column (toOriginalSpace) follows
    /// from the solve's own columns: scale * x[col] + constant when col >= 0,
    /// the constant when col == kImageConstant (presolve fixed it), and no
    /// such map when col == kImageNone (presolve substituted it through rows
    /// or bounds other than a doubleton equation, or split it with a
    /// duplicate column). Two top-level columns may share an image's column.
    static constexpr HighsInt kImageConstant = -1;
    static constexpr HighsInt kImageNone = -2;
    struct ColumnImage {
        HighsInt col;
        double scale;
        double constant;
    };

    /// The image of every top-level column in the solve's columns, through
    /// every sub-MIP level; identity at the top level. Reads the postsolve
    /// stacks as toOriginalSpace does.
    static std::vector<ColumnImage> columnImages(const HighsMipSolver& mipsolver);

    /// Guard against a reduction columnImages does not know (a HiGHS bump):
    /// clear the image of every top-level column whose image, evaluated at
    /// the solve's solution sol, differs from its expanded value
    /// (toOriginalSpace(sol)) by more than tol * (1 + |expanded value|).
    /// Returns the number of images cleared.
    static HighsInt clearMismatchedImages(std::vector<ColumnImage>& images,
                                          const std::vector<double>& sol,
                                          const std::vector<double>& expanded, double tol);

    /// Rewrite the cut idx . val <= rhs over top-level columns as one over
    /// the solve's columns, given the images columnImages returns: a fixed
    /// column moves to the rhs, and an affine one maps onto its column, with
    /// the coefficients of columns that map onto one column summed.
    /// Returns false, leaving the cut unspecified, when a cut column has no
    /// image: the cut cannot be stated exactly in the solve's columns.
    static bool toSolveSpace(const std::vector<ColumnImage>& images, std::vector<HighsInt>& idx,
                             std::vector<double>& val, double& rhs);

    /// Record that addIncumbent() rejected a solution because it violates the
    /// user lazy constraints. thread_local: the mark and its take happen in one
    /// synchronous call on one thread (HighsSearch::evaluateNode), on the main
    /// solver or on a tree worker.
    static void markIncumbentRejected();

    /// Return whether a rejection was recorded since the last call, and clear it.
    static bool takeIncumbentRejected();
};

#endif  // MIP_HIGHS_USER_SEPARATOR_H_
