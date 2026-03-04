#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <span>
#include <tuple>
#include <vector>

namespace cptp {

/// Simple directed graph with CSR storage for out-arcs and in-arcs.
/// Used by separators to build support graphs for max-flow.
class digraph {
 public:
    digraph() = default;

    int32_t num_nodes() const { return num_nodes_; }
    int32_t num_arcs() const { return num_arcs_; }

    std::span<const int32_t> out_arcs(int32_t v) const {
        return {out_arcs_.data() + out_begin_[v],
                out_arcs_.data() + out_begin_[v + 1]};
    }

    std::span<const int32_t> in_arcs(int32_t v) const {
        return {in_arcs_.data() + in_begin_[v],
                in_arcs_.data() + in_begin_[v + 1]};
    }

    int32_t arc_source(int32_t a) const { return sources_[a]; }
    int32_t arc_target(int32_t a) const { return targets_[a]; }

 private:
    friend class digraph_builder;
    int32_t num_nodes_ = 0;
    int32_t num_arcs_ = 0;
    std::vector<int32_t> sources_;
    std::vector<int32_t> targets_;
    std::vector<int32_t> out_begin_;  // CSR, size n+1
    std::vector<int32_t> out_arcs_;
    std::vector<int32_t> in_begin_;   // CSR, size n+1
    std::vector<int32_t> in_arcs_;
};

/// Builder for digraph. Collects arcs, then builds CSR.
class digraph_builder {
 public:
    explicit digraph_builder(int32_t n) : num_nodes_(n) {}

    void add_arc(int32_t u, int32_t v, double cap) {
        arc_sources_.push_back(u);
        arc_targets_.push_back(v);
        caps_.push_back(cap);
    }

    std::tuple<digraph, std::vector<double>> build() {
        digraph g;
        g.num_nodes_ = num_nodes_;
        g.num_arcs_ = static_cast<int32_t>(arc_sources_.size());
        g.sources_ = std::move(arc_sources_);
        g.targets_ = std::move(arc_targets_);

        // Sort arcs by source to build out-CSR
        std::vector<int32_t> perm(g.num_arcs_);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&](int32_t a, int32_t b) {
            return g.sources_[a] < g.sources_[b];
        });

        // Apply permutation
        auto permute = [&](auto& vec) {
            auto tmp = vec;
            for (int32_t i = 0; i < g.num_arcs_; ++i) vec[i] = tmp[perm[i]];
        };
        permute(g.sources_);
        permute(g.targets_);
        permute(caps_);

        // Build out-CSR
        g.out_begin_.assign(num_nodes_ + 1, 0);
        for (int32_t a = 0; a < g.num_arcs_; ++a) {
            g.out_begin_[g.sources_[a] + 1]++;
        }
        for (int32_t i = 0; i < num_nodes_; ++i) {
            g.out_begin_[i + 1] += g.out_begin_[i];
        }
        g.out_arcs_.resize(g.num_arcs_);
        std::iota(g.out_arcs_.begin(), g.out_arcs_.end(), 0);

        // Build in-CSR
        g.in_begin_.assign(num_nodes_ + 1, 0);
        for (int32_t a = 0; a < g.num_arcs_; ++a) {
            g.in_begin_[g.targets_[a] + 1]++;
        }
        for (int32_t i = 0; i < num_nodes_; ++i) {
            g.in_begin_[i + 1] += g.in_begin_[i];
        }
        g.in_arcs_.resize(g.num_arcs_);
        std::vector<int32_t> pos(num_nodes_, 0);
        for (int32_t a = 0; a < g.num_arcs_; ++a) {
            int32_t t = g.targets_[a];
            g.in_arcs_[g.in_begin_[t] + pos[t]++] = a;
        }

        return {std::move(g), std::move(caps_)};
    }

 private:
    int32_t num_nodes_;
    std::vector<int32_t> arc_sources_;
    std::vector<int32_t> arc_targets_;
    std::vector<double> caps_;
};

}  // namespace cptp
