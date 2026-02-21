#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include "core/digraph.h"

namespace rcspp {

/// Dinitz max-flow on a digraph with double capacities.
class dinitz {
 public:
    dinitz(const digraph& g, const std::vector<double>& cap, int32_t s, int32_t t)
        : g_(g), cap_(cap), source_(s), target_(t),
          flow_(g.num_arcs(), 0.0),
          level_(g.num_nodes()),
          iter_(g.num_nodes()) {
        // Build reverse-arc mapping: for arc a, rev_[a] is its reverse arc.
        // We assume arcs come in forward/reverse pairs OR we need to find them.
        // Actually, our support graphs are built with explicit forward+reverse arcs.
        // We need to find the reverse of each arc efficiently.
        int32_t m = g.num_arcs();
        rev_.assign(m, -1);

        // For each arc (u,v), find the arc (v,u) as its reverse
        // Build a map: (u,v) -> arc index
        for (int32_t a = 0; a < m; ++a) {
            if (rev_[a] >= 0) continue;
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
        std::fill(flow_.begin(), flow_.end(), 0.0);

        while (bfs()) {
            std::fill(iter_.begin(), iter_.end(), 0);
            double f;
            while ((f = dfs(source_, std::numeric_limits<double>::infinity())) > 0) {
                flow_value_ += f;
            }
        }
        return *this;
    }

    double flow_value() const { return flow_value_; }

    /// Returns true if v is on the source side of the min-cut.
    bool on_source_side(int32_t v) const {
        return level_[v] >= 0;
    }

 private:
    bool bfs() {
        std::fill(level_.begin(), level_.end(), -1);
        std::queue<int32_t> q;
        level_[source_] = 0;
        q.push(source_);

        while (!q.empty()) {
            int32_t v = q.front(); q.pop();
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

    double dfs(int32_t v, double pushed) {
        if (v == target_) return pushed;

        // Try forward arcs
        auto out = g_.out_arcs(v);
        auto in = g_.in_arcs(v);
        int32_t total = static_cast<int32_t>(out.size() + in.size());

        for (int32_t& i = iter_[v]; i < total; ++i) {
            if (i < static_cast<int32_t>(out.size())) {
                int32_t a = out[i];
                int32_t w = g_.arc_target(a);
                double r = residual(a);
                if (r > 1e-12 && level_[w] == level_[v] + 1) {
                    double f = dfs(w, std::min(pushed, r));
                    if (f > 1e-12) {
                        flow_[a] += f;
                        return f;
                    }
                }
            } else {
                // Reverse direction: push back flow on in-arc
                int32_t idx = i - static_cast<int32_t>(out.size());
                int32_t a = in[idx];
                int32_t w = g_.arc_source(a);
                if (flow_[a] > 1e-12 && level_[w] == level_[v] + 1) {
                    double f = dfs(w, std::min(pushed, flow_[a]));
                    if (f > 1e-12) {
                        flow_[a] -= f;
                        return f;
                    }
                }
            }
        }

        return 0.0;
    }

    double residual(int32_t a) const {
        return cap_[a] - flow_[a];
    }

    const digraph& g_;
    const std::vector<double>& cap_;
    int32_t source_, target_;
    std::vector<double> flow_;
    std::vector<int32_t> level_;
    std::vector<int32_t> iter_;
    std::vector<int32_t> rev_;
    double flow_value_ = 0.0;
};

}  // namespace rcspp
