#include "mip/HighsUserSeparator.h"
#include "mip/HighsMipSolverData.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

// Static member definitions — single copy linked into libhighs.
HighsUserSeparator::Callback HighsUserSeparator::callback_;
HighsUserSeparator::FeasibilityCheck HighsUserSeparator::feasibility_check_;
HighsUserSeparator::BranchingCallback HighsUserSeparator::branching_callback_;
HighsUserSeparator::StrongBranchConfig HighsUserSeparator::strong_branch_config_;
std::mutex HighsUserSeparator::mutex_;
std::unordered_map<std::string, HighsUserSeparator::HyperplanePseudocost>
    HighsUserSeparator::hyperplane_pseudocosts_;
std::vector<HighsInt> HighsUserSeparator::hp_stack_;
std::vector<HighsUserSeparator::HyperplaneData> HighsUserSeparator::hp_data_pool_;

void HighsUserSeparator::separateLpSolution(
    HighsLpRelaxation& lpRelaxation,
    HighsLpAggregator& /*lpAggregator*/,
    HighsTransformedLp& /*transLp*/,
    HighsCutPool& cutpool) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (callback_) {
    callback_(lpRelaxation, cutpool, lpRelaxation.getMipSolver());
  }
}

void HighsUserSeparator::setCallback(Callback cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = std::move(cb);
}

void HighsUserSeparator::setFeasibilityCheck(FeasibilityCheck cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  feasibility_check_ = std::move(cb);
}

void HighsUserSeparator::clearCallback() {
  std::lock_guard<std::mutex> lock(mutex_);
  callback_ = nullptr;
  feasibility_check_ = nullptr;
  branching_callback_ = nullptr;
  strong_branch_config_ = StrongBranchConfig{};
  hyperplane_pseudocosts_.clear();
  clearStack();
}

void HighsUserSeparator::setBranchingCallback(BranchingCallback cb) {
  std::lock_guard<std::mutex> lock(mutex_);
  branching_callback_ = std::move(cb);
}

const HighsUserSeparator::BranchingCallback&
HighsUserSeparator::getBranchingCallback() {
  return branching_callback_;
}

void HighsUserSeparator::setStrongBranchConfig(StrongBranchConfig cfg) {
  strong_branch_config_ = cfg;
}

const HighsUserSeparator::StrongBranchConfig&
HighsUserSeparator::getStrongBranchConfig() {
  return strong_branch_config_;
}

bool HighsUserSeparator::hasCallback() {
  std::lock_guard<std::mutex> lock(mutex_);
  return static_cast<bool>(callback_);
}

bool HighsUserSeparator::isFeasible(const std::vector<double>& sol) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!feasibility_check_) return true;
  return feasibility_check_(sol);
}

// ── Hyperplane branching stack/pool management ──

HighsInt HighsUserSeparator::allocHyperplaneData() {
  HighsInt idx = static_cast<HighsInt>(hp_data_pool_.size());
  hp_data_pool_.emplace_back();
  return idx;
}

HighsUserSeparator::HyperplaneData&
HighsUserSeparator::getHyperplaneData(HighsInt idx) {
  return hp_data_pool_[idx];
}

HighsInt HighsUserSeparator::poolSize() {
  return static_cast<HighsInt>(hp_data_pool_.size());
}

HighsInt HighsUserSeparator::stackSize() {
  return static_cast<HighsInt>(hp_stack_.size());
}

void HighsUserSeparator::pushBranching(HighsLpRelaxation& lp,
                                        const HyperplaneData& data,
                                        double lo, double hi) {
  HighsInt nnz = static_cast<HighsInt>(data.indices.size());
  lp.addBranchingRow(nnz, data.indices.data(), data.values.data(), lo, hi);
  hp_stack_.push_back(1);
}

void HighsUserSeparator::restoreBranching(HighsLpRelaxation& lp,
                                           HighsInt target_pos) {
  while (static_cast<HighsInt>(hp_stack_.size()) > target_pos) {
    lp.removeLastBranchingRow();
    hp_stack_.pop_back();
  }
}

void HighsUserSeparator::clearStack() {
  hp_stack_.clear();
  hp_data_pool_.clear();
}

// ── Candidate evaluation ──

int HighsUserSeparator::evaluateCandidates(
    std::vector<HyperplaneCandidate>& candidates,
    double var_score,
    const HighsMipSolver& mipsolver,
    int depth) {
  const auto& cfg = strong_branch_config_;
  const auto& sol = mipsolver.mipdata_->lp.getSolution().col_value;

  // Compute LHS and fractionality for each candidate, filter fractional ones
  struct CandInfo {
    int idx;
    double lhs;
    double frac;       // fractional part
    double dual_score;
  };
  std::vector<CandInfo> fractional;

  // Get reduced costs for dual pre-scoring
  const auto& col_dual = mipsolver.mipdata_->lp.getLpSolver()
                              .getSolution().col_dual;

  for (int c = 0; c < (int)candidates.size(); ++c) {
    double lhs = 0.0;
    for (int j = 0; j < (int)candidates[c].indices.size(); ++j)
      lhs += candidates[c].values[j] * sol[candidates[c].indices[j]];
    double frac = lhs - std::floor(lhs);
    if (frac < 0.01 || frac > 0.99) continue;

    // Tier 1: dual pre-score
    double ds = 0.0;
    for (int j = 0; j < (int)candidates[c].indices.size(); ++j) {
      ds += std::abs(candidates[c].values[j] *
                     col_dual[candidates[c].indices[j]]);
    }
    ds *= std::min(frac, 1.0 - frac);  // weight by fractionality

    fractional.push_back({c, lhs, frac, ds});
  }

  if (fractional.empty()) return -1;

  // Sort by dual score descending, keep top-k
  std::sort(fractional.begin(), fractional.end(),
            [](const CandInfo& a, const CandInfo& b) {
              return a.dual_score > b.dual_score;
            });
  if ((int)fractional.size() > cfg.max_sb_candidates)
    fractional.resize(cfg.max_sb_candidates);

  // Tier 2 & 3: pseudocost lookup or sequential SB
  double best_score = -1.0;
  int best_idx = -1;

  bool do_sb = (cfg.max_depth > 0 && depth <= cfg.max_depth);
  auto& lp = mipsolver.mipdata_->lp;
  auto& lpsolver = lp.getLpSolver();
  double parentObj = lp.getObjective();

  // Save settings for SB and force simplex
  HighsInt orig_iter_limit;
  lpsolver.getOptionValue("simplex_iteration_limit", orig_iter_limit);
  std::string orig_solver;
  lpsolver.getOptionValue("solver", orig_solver);
  lpsolver.setOptionValue("solver", "simplex");

  for (auto& ci : fractional) {
    auto& cand = candidates[ci.idx];
    auto& pc = hyperplane_pseudocosts_[cand.name];

    // Tier 2: pseudocost reliable → use estimate
    if (pc.reliable(cfg.min_reliable)) {
      double score = pc.score(ci.lhs);
      if (score > best_score) {
        best_score = score;
        best_idx = ci.idx;
      }
      continue;
    }

    // Tier 3: sequential SB (only within depth limit)
    if (!do_sb) {
      // Beyond depth limit with unreliable pseudocost: use dual score as proxy
      double score = ci.dual_score;
      if (score > best_score) {
        best_score = score;
        best_idx = ci.idx;
      }
      continue;
    }

    lpsolver.setOptionValue("simplex_iteration_limit", cfg.iter_limit);

    // SB via dynamic branching rows
    HyperplaneData tmp_data;
    tmp_data.indices = cand.indices;
    tmp_data.values = cand.values;

    // Down branch: lhs <= floor(lhs)
    lp.addBranchingRow(
        static_cast<HighsInt>(tmp_data.indices.size()),
        tmp_data.indices.data(), tmp_data.values.data(),
        -kHighsInf, std::floor(ci.lhs));
    lpsolver.run();
    auto status_down = lpsolver.getModelStatus();
    double down_lb = (status_down == HighsModelStatus::kInfeasible)
        ? kHighsInf : lpsolver.getObjectiveValue();
    lp.removeLastBranchingRow();

    // Up branch: lhs >= ceil(lhs)
    lp.addBranchingRow(
        static_cast<HighsInt>(tmp_data.indices.size()),
        tmp_data.indices.data(), tmp_data.values.data(),
        std::ceil(ci.lhs), kHighsInf);
    lpsolver.run();
    auto status_up = lpsolver.getModelStatus();
    double up_lb = (status_up == HighsModelStatus::kInfeasible)
        ? kHighsInf : lpsolver.getObjectiveValue();
    lp.removeLastBranchingRow();

    // Update pseudocost
    double di = std::max(0.0, down_lb - parentObj);
    double ui = std::max(0.0, up_lb - parentObj);
    pc.addObservation(ci.lhs, di, false);
    pc.addObservation(ci.lhs, ui, true);

    double score = std::min(di, ui) * std::max(di, ui);
    if (score > best_score) {
      best_score = score;
      best_idx = ci.idx;
    }
  }

  // Restore settings and re-solve to recover original LP solution state
  lpsolver.setOptionValue("simplex_iteration_limit", orig_iter_limit);
  lpsolver.setOptionValue("solver", orig_solver);
  if (do_sb) {
    // SB trials may have modified the LP basis/solution — re-solve to recover
    lpsolver.run();
  }

  // Compare best hyperplane score vs variable branching score
  if (best_idx >= 0 && best_score > var_score) {
    return best_idx;
  }
  return -1;
}
