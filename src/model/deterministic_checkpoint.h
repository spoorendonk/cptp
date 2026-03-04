#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <utility>
#include <vector>

namespace cptp::model_detail {

struct DssrEpochUpdate {
    int32_t epoch = 0;
    int32_t ng_size = 1;
    std::vector<double> fwd;
    std::vector<double> bwd;
    bool elementary_path_found = false;
    double elementary_path_cost = std::numeric_limits<double>::infinity();
    std::vector<int32_t> elementary_path;
};

class DeterministicCheckpointClock {
 public:
    uint64_t checkpoint() {
        std::lock_guard<std::mutex> lock(mutex_);
        return ++value_;
    }

    uint64_t value() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return value_;
    }

 private:
    mutable std::mutex mutex_;
    uint64_t value_ = 0;
};

class DssrAutoStopTracker {
 public:
    DssrAutoStopTracker(int32_t min_epochs_before_stop,
                        int32_t no_progress_epoch_limit)
        : min_epochs_before_stop_(std::max<int32_t>(1, min_epochs_before_stop)),
          no_progress_epoch_limit_(std::max<int32_t>(1, no_progress_epoch_limit)) {}

    bool observe_stage(bool improved, bool checkpoint_active, int32_t produced_epochs) {
        had_checkpoint_activity_ = had_checkpoint_activity_ || checkpoint_active;
        if (improved) {
            no_progress_epochs_ = 0;
        } else {
            ++no_progress_epochs_;
        }
        if (produced_epochs >= min_epochs_before_stop_ && !had_checkpoint_activity_) {
            return true;
        }
        if (no_progress_epochs_ >= no_progress_epoch_limit_) {
            return true;
        }
        return false;
    }

    int32_t no_progress_epochs() const { return no_progress_epochs_; }
    bool had_checkpoint_activity() const { return had_checkpoint_activity_; }

 private:
    int32_t min_epochs_before_stop_ = 1;
    int32_t no_progress_epoch_limit_ = 1;
    int32_t no_progress_epochs_ = 0;
    bool had_checkpoint_activity_ = false;
};

// Thread-safe queue keyed by epoch number. Commits are performed in epoch order.
class DssrEpochQueue {
 public:
    explicit DssrEpochQueue(int32_t first_epoch)
        : next_epoch_(first_epoch) {}

    bool enqueue(DssrEpochUpdate update) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto [it, inserted] = pending_.emplace(update.epoch, std::move(update));
        if (!inserted) {
            it->second = std::move(update);
        }
        cv_.notify_all();
        return inserted;
    }

    void mark_producer_done() {
        std::lock_guard<std::mutex> lock(mutex_);
        producer_done_ = true;
        cv_.notify_all();
    }

    template <typename CommitFn>
    bool commit_next_blocking(CommitFn&& commit) {
        std::optional<DssrEpochUpdate> next;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            cv_.wait(lock, [&] {
                return producer_done_ || pending_.contains(next_epoch_);
            });
            auto it = pending_.find(next_epoch_);
            if (it == pending_.end()) return false;
            next = std::move(it->second);
            pending_.erase(it);
            ++next_epoch_;
        }
        commit(std::move(*next));
        return true;
    }

    template <typename CommitFn>
    int64_t commit_ready(CommitFn&& commit) {
        int64_t committed = 0;
        while (true) {
            auto next = take_ready_locked();
            if (!next.has_value()) break;
            commit(std::move(*next));
            committed++;
        }
        return committed;
    }

    template <typename CommitFn>
    int64_t commit_all_ordered(CommitFn&& commit) {
        int64_t committed = 0;
        while (true) {
            auto next = take_smallest_locked();
            if (!next.has_value()) break;
            commit(std::move(*next));
            committed++;
        }
        return committed;
    }

 private:
    std::optional<DssrEpochUpdate> take_ready_locked() {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = pending_.find(next_epoch_);
        if (it == pending_.end()) return std::nullopt;
        DssrEpochUpdate out = std::move(it->second);
        pending_.erase(it);
        ++next_epoch_;
        return out;
    }

    std::optional<DssrEpochUpdate> take_smallest_locked() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (pending_.empty()) return std::nullopt;
        auto it = pending_.begin();
        DssrEpochUpdate out = std::move(it->second);
        next_epoch_ = std::max(next_epoch_, out.epoch + 1);
        pending_.erase(it);
        return out;
    }

    std::mutex mutex_;
    std::condition_variable cv_;
    std::map<int32_t, DssrEpochUpdate> pending_;
    int32_t next_epoch_ = 1;
    bool producer_done_ = false;
};

}  // namespace cptp::model_detail
