#ifndef MIP_HIGHS_USER_HOOKS_H_
#define MIP_HIGHS_USER_HOOKS_H_

#include <functional>
#include <memory>
#include <vector>

class HighsCutPool;
class HighsDomain;
class HighsLpRelaxation;
class HighsMipSolver;

/// The user callbacks one MIP solve runs with: the separator, the lazy
/// feasibility check, the domain propagator and the root task observer.
///
/// The static setters of HighsUserSeparator and HighsUserPropagator fill a
/// thread_local pending set on the thread that calls Highs::run(), so
/// concurrent Highs instances on different threads stay apart. A top-level
/// HighsMipSolver::run() captures that set into HighsMipSolver::user_hooks_
/// before any task is spawned, and a sub-MIP copies its parent's: every
/// worker thread of HiGHS's parallel tree search, and every sub-MIP started
/// on one, then runs the same hooks.
struct HighsUserHooks {
    using Separator = std::function<void(const HighsLpRelaxation& lpRelaxation,
                                         HighsCutPool& cutpool, const HighsMipSolver& mipsolver)>;
    /// Returns true if the solution violates no lazy constraint.
    using FeasibilityCheck = std::function<bool(const std::vector<double>& sol)>;
    using Propagator = std::function<void(HighsDomain& domain, const HighsMipSolver& mipsolver,
                                          const HighsLpRelaxation& lp)>;

    Separator separator;
    FeasibilityCheck feasibility_check;
    Propagator propagator;

    /// Limits of the main MIP's root cut loop (evaluateRootNode): its round
    /// cap and its stall limit (rounds without progress before it stops).
    /// kRootLimitHighs keeps HiGHS's own value, kRootLimitUnlimited lifts the
    /// limit, a positive value replaces HiGHS's. A sub-MIP keeps HiGHS's.
    static constexpr int kRootLimitHighs = -1;
    static constexpr int kRootLimitUnlimited = 0;
    int root_sepa_rounds = kRootLimitHighs;
    int root_sepa_stall = kRootLimitHighs;

    /// Told +1 when a root task HiGHS spawns beside its root node (symmetry
    /// detection, the analytic centre), of the solve or one of its sub-MIPs,
    /// starts on a thread other than the one that spawned it, and -1 when
    /// that task ends: the threads those tasks keep busy while the callbacks
    /// run. Called on the task's thread. Null leaves them uncounted.
    using RootTaskObserver = std::function<void(int delta)>;
    RootTaskObserver root_task_observer;

    /// The pending set of the calling thread (defined in HighsUserSeparator.cpp).
    static HighsUserHooks& pending();

    /// A copy of the calling thread's pending set, or null when it has no
    /// callback and HiGHS's root limits.
    static std::shared_ptr<const HighsUserHooks> capture();
};

#endif  // MIP_HIGHS_USER_HOOKS_H_
