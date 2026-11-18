#pragma once

#include "core/digraph.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

namespace cptp {

/// Dinitz max-flow on a digraph with double capacities.
class dinitz {
public:
    dinitz(const digraph& g, const std::vector<double>& cap, int32_t s, int32_t t)
        : g_(g),
          cap_(cap),
          source_(s),
          target_(t),
          flow_(g.num_arcs(), 0.0),
          level_(g.num_nodes()),
          iter_(g.num_nodes()) {
        // Build reverse-arc mapping: for arc a, rev_[a] is its reverse arc.
        // We assume arcs come in forward/reverse pairs OR we need to find them.
        // Actually, our support graphs are built with explicit forward+reverse
        // arcs. We need to find the reverse of each arc efficiently.
        int32_t m = g.num_arcs();
        rev_.assign(m, -1);

        // For each arc (u,v), find the arc (v,u) as its reverse
        // Build a map: (u,v) -> arc index
        for (int32_t a = 0; a < m; ++a) {
            if (rev_[a] >= 0) {
                continue;
            }
            int32_t u = g.arc_source(a);
            int32_t v = g.arc_target(a);
            // Find reverse arc (v,u)
            for (int32_t b : g.out_arcs(v)) {
                if (g.arc_target(b) == u && rev_[b] < 0) {
                    rev_[a] = b;
                    rev_[b] = a;
                    break;
                }
            }
        }
    }

    dinitz& run() {
        flow_value_ = 0.0;
        std::ranges::fill(flow_, 0.0);

        while (bfs()) {
            std::ranges::fill(iter_, 0);
            double f;
            while ((f = dfs()) > 0) {
                flow_value_ += f;
            }
        }
        return *this;
    }

    [[nodiscard]] double flow_value() const { return flow_value_; }

    /// Returns true if v is on the source side of the min-cut.
    [[nodiscard]] bool on_source_side(int32_t v) const { return level_[v] >= 0; }

private:
    bool bfs() {
        std::ranges::fill(level_, -1);
        std::queue<int32_t> q;
        level_[source_] = 0;
        q.push(source_);

        while (!q.empty()) {
            int32_t v = q.front();
            q.pop();
            for (int32_t a : g_.out_arcs(v)) {
                int32_t w = g_.arc_target(a);
                if (level_[w] < 0 && residual(a) > 1e-12) {
                    level_[w] = level_[v] + 1;
                    q.push(w);
                }
            }
            // Also check reverse arcs (flow can be pushed back)
            for (int32_t a : g_.in_arcs(v)) {
                int32_t w = g_.arc_source(a);
                if (level_[w] < 0 && flow_[a] > 1e-12) {
                    level_[w] = level_[v] + 1;
                    q.push(w);
                }
            }
        }

        return level_[target_] >= 0;
    }

    /// One augmenting path in the level graph, found with the current-arc
    /// pointers `iter_`; returns its bottleneck, or 0 when none is left.
    ///
    /// Iterative, with the path on the heap: the recursion it replaces went
    /// one frame per BFS level, up to n-1 deep on a path-like support, and
    /// overflowed the 512 KB stack of a macOS secondary thread (HiGHS worker
    /// threads run this) at a few thousand nodes (#161). It finds the same
    /// path the recursion did: an arc is skipped, advancing `iter_[v]`, only
    /// when the search below it is exhausted, and the arc that led to the
    /// target stays current.
    double dfs() {
        path_.clear();
        int32_t v = source_;
        while (v != target_) {
            auto out = g_.out_arcs(v);
            auto in = g_.in_arcs(v);
            auto num_out = static_cast<int32_t>(out.size());
            auto total = static_cast<int32_t>(out.size() + in.size());

            int32_t next = -1;
            for (int32_t& i = iter_[v]; i < total; ++i) {
                if (i < num_out) {
                    int32_t a = out[i];
                    int32_t w = g_.arc_target(a);
                    if (residual(a) > 1e-12 && level_[w] == level_[v] + 1) {
                        path_.push_back({.from = v, .arc = a, .forward = true});
                        next = w;
                        break;
                    }
                } else {
                    // Reverse direction: push back flow on in-arc
                    int32_t a = in[i - num_out];
                    int32_t w = g_.arc_source(a);
                    if (flow_[a] > 1e-12 && level_[w] == level_[v] + 1) {
                        path_.push_back({.from = v, .arc = a, .forward = false});
                        next = w;
                        break;
                    }
                }
            }

            if (next >= 0) {
                v = next;
                continue;
            }
            // v is a dead end: retreat and skip the arc that led here.
            if (path_.empty()) {
                return 0.0;
            }
            v = path_.back().from;
            path_.pop_back();
            ++iter_[v];
        }

        double f = std::numeric_limits<double>::infinity();
        for (const Step& s : path_) {
            f = std::min(f, s.forward ? residual(s.arc) : flow_[s.arc]);
        }
        for (const Step& s : path_) {
            flow_[s.arc] += s.forward ? f : -f;
        }
        return f;
    }

    [[nodiscard]] double residual(int32_t a) const { return cap_[a] - flow_[a]; }

    /// An arc of the current augmenting path, taken out of node `from`:
    /// forward along its residual, or backward against its flow.
    struct Step {
        int32_t from;
        int32_t arc;
        bool forward;
    };

    const digraph& g_;
    const std::vector<double>& cap_;
    int32_t source_, target_;
    std::vector<double> flow_;
    std::vector<int32_t> level_;
    std::vector<int32_t> iter_;
    std::vector<int32_t> rev_;
    std::vector<Step> path_;
    double flow_value_ = 0.0;
};

}  // namespace cptp
