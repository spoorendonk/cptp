#include "mip/HighsUserSeparator.h"

#include "mip/HighsMipSolverData.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace {
// The calling thread's pending hooks, captured by the next top-level
// HighsMipSolver::run() on this thread.
thread_local HighsUserHooks pending_hooks_;
thread_local bool incumbent_rejected_ = false;
}  // namespace

HighsUserHooks& HighsUserHooks::pending() {
    return pending_hooks_;
}

std::shared_ptr<const HighsUserHooks> HighsUserHooks::capture() {
    const HighsUserHooks& hooks = pending_hooks_;
    if (!hooks.separator && !hooks.feasibility_check && !hooks.propagator &&
        !hooks.root_task_observer && hooks.root_sepa_rounds == HighsUserHooks::kRootLimitHighs &&
        hooks.root_sepa_stall == HighsUserHooks::kRootLimitHighs)
        return nullptr;
    return std::make_shared<const HighsUserHooks>(hooks);
}

void HighsUserSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& /*lpAggregator*/,
                                            HighsTransformedLp& /*transLp*/,
                                            HighsCutPool& cutpool) {
    const HighsMipSolver& mipsolver = lpRelaxation.getMipSolver();
    if (hasCallback(mipsolver)) {
        mipsolver.user_hooks_->separator(lpRelaxation, cutpool, mipsolver);
    }
}

void HighsUserSeparator::setCallback(Callback cb) {
    pending_hooks_.separator = std::move(cb);
}

void HighsUserSeparator::setFeasibilityCheck(FeasibilityCheck cb) {
    pending_hooks_.feasibility_check = std::move(cb);
}

void HighsUserSeparator::setRootSeparationLimits(int max_rounds, int stall_limit) {
    pending_hooks_.root_sepa_rounds = max_rounds;
    pending_hooks_.root_sepa_stall = stall_limit;
}

void HighsUserSeparator::setRootTaskObserver(HighsUserHooks::RootTaskObserver observer) {
    pending_hooks_.root_task_observer = std::move(observer);
}

void HighsUserSeparator::clearCallback() {
    pending_hooks_.separator = nullptr;
    pending_hooks_.feasibility_check = nullptr;
    pending_hooks_.root_task_observer = nullptr;
    pending_hooks_.root_sepa_rounds = HighsUserHooks::kRootLimitHighs;
    pending_hooks_.root_sepa_stall = HighsUserHooks::kRootLimitHighs;
}

namespace {
// The solve's root task observer, or null.
const HighsUserHooks::RootTaskObserver* rootTaskObserver(const HighsMipSolver& mipsolver) {
    if (!mipsolver.user_hooks_ || !mipsolver.user_hooks_->root_task_observer)
        return nullptr;
    return &mipsolver.user_hooks_->root_task_observer;
}
}  // namespace

void HighsUserSeparator::rootTaskSpawned(const HighsMipSolver& mipsolver) {
    const auto* observer = rootTaskObserver(mipsolver);
    if (!observer)
        return;
    mipsolver.root_tasks_pending_.fetch_add(1, std::memory_order_relaxed);
    (*observer)(+1, 0);
}

void HighsUserSeparator::rootTasksJoined(const HighsMipSolver& mipsolver) {
    const auto* observer = rootTaskObserver(mipsolver);
    if (!observer)
        return;
    // Joined: every task that ran has ended and released itself.
    const int unrun = mipsolver.root_tasks_pending_.exchange(0, std::memory_order_relaxed);
    if (unrun != 0)
        (*observer)(-unrun, 0);
}

HighsUserSeparator::RootTaskScope::RootTaskScope(const HighsMipSolver& mipsolver,
                                                 std::thread::id owner)
    : mipsolver_(mipsolver) {
    const auto* observer = rootTaskObserver(mipsolver_);
    if (!observer || std::this_thread::get_id() == owner)
        return;
    running_ = true;
    (*observer)(0, +1);
}

HighsUserSeparator::RootTaskScope::~RootTaskScope() {
    const auto* observer = rootTaskObserver(mipsolver_);
    if (!observer)
        return;
    mipsolver_.root_tasks_pending_.fetch_sub(1, std::memory_order_relaxed);
    (*observer)(-1, running_ ? -1 : 0);
}

namespace {
// A root limit of the solve's hooks in place of HiGHS's own (HighsUserHooks).
HighsInt rootLimit(const HighsMipSolver& mipsolver, int HighsUserHooks::* limit,
                   HighsInt highs_value) {
    if (mipsolver.submip || !mipsolver.user_hooks_)
        return highs_value;
    const int value = (*mipsolver.user_hooks_).*limit;
    if (value == HighsUserHooks::kRootLimitUnlimited)
        return kHighsIInf;
    if (value > 0)
        return value;
    return highs_value;
}
}  // namespace

HighsInt HighsUserSeparator::rootSepaRounds(const HighsMipSolver& mipsolver,
                                            HighsInt highs_rounds) {
    return rootLimit(mipsolver, &HighsUserHooks::root_sepa_rounds, highs_rounds);
}

HighsInt HighsUserSeparator::rootSepaStall(const HighsMipSolver& mipsolver, HighsInt highs_stall) {
    return rootLimit(mipsolver, &HighsUserHooks::root_sepa_stall, highs_stall);
}

bool HighsUserSeparator::hasCallback() {
    return static_cast<bool>(pending_hooks_.separator);
}

bool HighsUserSeparator::isFeasible(const std::vector<double>& sol) {
    if (!pending_hooks_.feasibility_check)
        return true;
    return pending_hooks_.feasibility_check(sol);
}

bool HighsUserSeparator::hasCallback(const HighsMipSolver& mipsolver) {
    return mipsolver.user_hooks_ && mipsolver.user_hooks_->separator;
}

bool HighsUserSeparator::isFeasible(const HighsMipSolver& mipsolver,
                                    const std::vector<double>& sol) {
    if (!mipsolver.user_hooks_ || !mipsolver.user_hooks_->feasibility_check)
        return true;
    const auto& check = mipsolver.user_hooks_->feasibility_check;
    if (!mipsolver.submip)
        return check(sol);
    return check(toOriginalSpace(mipsolver, sol));
}

std::vector<double> HighsUserSeparator::toOriginalSpace(const HighsMipSolver& mipsolver,
                                                        const std::vector<double>& sol) {
    std::vector<double> col_value = sol;
    // Each sub-MIP's model is its parent's presolved model, so one undo per
    // level. The walk stops below the top-level solve, whose presolve is off;
    // a sub-MIP runs on one thread, inside its parent's call, so the stacks
    // read here are not in use elsewhere.
    for (const HighsMipSolver* level = &mipsolver; level->submip; level = level->parent_mip_) {
        // A fresh solution per level: undo sizes row_value for one level only.
        HighsSolution expanded;
        expanded.col_value = std::move(col_value);
        expanded.value_valid = true;
        level->mipdata_->postSolveStack.undoPrimal(*level->options_mip_, expanded);
        col_value = std::move(expanded.col_value);
    }
    return col_value;
}

std::vector<HighsInt> HighsUserSeparator::originalColumns(const HighsMipSolver& mipsolver,
                                                          HighsInt num_col) {
    std::vector<HighsInt> cols(static_cast<size_t>(num_col));
    std::iota(cols.begin(), cols.end(), HighsInt{0});
    for (const HighsMipSolver* level = &mipsolver; level->submip; level = level->parent_mip_) {
        const auto& stack = level->mipdata_->postSolveStack;
        for (HighsInt& col : cols)
            col = stack.getOrigColIndex(col);
    }
    return cols;
}

std::vector<HighsUserSeparator::ColumnImage> HighsUserSeparator::columnImages(
    const HighsMipSolver& mipsolver) {
    using Stack = presolve::HighsPostsolveStack;
    std::vector<const HighsMipSolver*> levels;
    for (const HighsMipSolver* level = &mipsolver; level->submip; level = level->parent_mip_)
        levels.push_back(level);

    // Start with the identity on the top-level model's columns: the input of
    // the outermost sub-MIP level, or the solve's own at the top level.
    const HighsInt num_top_col = levels.empty()
                                     ? mipsolver.model_->num_col_
                                     : levels.back()->mipdata_->postSolveStack.getOrigNumCol();
    std::vector<ColumnImage> images(static_cast<size_t>(num_top_col));
    for (HighsInt col = 0; col < num_top_col; ++col)
        images[static_cast<size_t>(col)] = {col, 1.0, 0.0};

    // Outermost level first: each maps its input columns, the columns the
    // images so far land on, into its presolved ones.
    for (auto it = levels.rbegin(); it != levels.rend(); ++it) {
        const std::vector<Stack::CptpColumnImage> level_images =
            (*it)->mipdata_->postSolveStack.cptpColumnImages();
        for (ColumnImage& image : images) {
            if (image.col < 0)
                continue;
            const Stack::CptpColumnImage& inner = level_images[static_cast<size_t>(image.col)];
            if (inner.col == Stack::kCptpImageNone) {
                image = {kImageNone, 0.0, 0.0};
            } else if (inner.col == Stack::kCptpImageConstant) {
                image = {kImageConstant, 0.0, image.constant + image.scale * inner.constant};
            } else {
                image = {inner.col, image.scale * inner.scale,
                         image.constant + image.scale * inner.constant};
            }
        }
    }
    return images;
}

HighsInt HighsUserSeparator::clearMismatchedImages(std::vector<ColumnImage>& images,
                                                   const std::vector<double>& sol,
                                                   const std::vector<double>& expanded,
                                                   double tol) {
    HighsInt cleared = 0;
    for (size_t j = 0; j < images.size(); ++j) {
        ColumnImage& image = images[j];
        if (image.col == kImageNone)
            continue;
        bool match = j < expanded.size();
        if (match) {
            double value = image.constant;
            if (image.col != kImageConstant) {
                match = static_cast<size_t>(image.col) < sol.size();
                if (match)
                    value += image.scale * sol[static_cast<size_t>(image.col)];
            }
            match = match && std::abs(value - expanded[j]) <= tol * (1.0 + std::abs(expanded[j]));
        }
        if (!match) {
            image = {kImageNone, 0.0, 0.0};
            ++cleared;
        }
    }
    return cleared;
}

bool HighsUserSeparator::toSolveSpace(const std::vector<ColumnImage>& images,
                                      std::vector<HighsInt>& idx, std::vector<double>& val,
                                      double& rhs) {
    std::vector<std::pair<HighsInt, double>> terms;
    terms.reserve(idx.size());
    for (size_t k = 0; k < idx.size(); ++k) {
        const HighsInt col = idx[k];
        if (col < 0 || static_cast<size_t>(col) >= images.size())
            return false;
        const ColumnImage& image = images[static_cast<size_t>(col)];
        if (image.col == kImageNone)
            return false;
        rhs -= val[k] * image.constant;
        if (image.col != kImageConstant)
            terms.emplace_back(image.col, val[k] * image.scale);
    }
    // Two top-level columns may land on one column (a doubleton equation's
    // pair): sum their coefficients, and drop a term that cancels.
    // HiGHS builds as C++11: no ranges or structured bindings here.
    std::stable_sort(terms.begin(), terms.end(),
                     [](const std::pair<HighsInt, double>& a,
                        const std::pair<HighsInt, double>& b) { return a.first < b.first; });
    idx.clear();
    val.clear();
    for (const std::pair<HighsInt, double>& term : terms) {
        if (!idx.empty() && idx.back() == term.first) {
            val.back() += term.second;
        } else {
            idx.push_back(term.first);
            val.push_back(term.second);
        }
    }
    size_t kept = 0;
    for (size_t k = 0; k < idx.size(); ++k) {
        if (val[k] == 0.0)
            continue;
        idx[kept] = idx[k];
        val[kept] = val[k];
        ++kept;
    }
    idx.resize(kept);
    val.resize(kept);
    return true;
}

void HighsUserSeparator::markIncumbentRejected() {
    incumbent_rejected_ = true;
}

bool HighsUserSeparator::takeIncumbentRejected() {
    const bool rejected = incumbent_rejected_;
    incumbent_rejected_ = false;
    return rejected;
}
