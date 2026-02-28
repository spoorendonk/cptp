#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "core/problem.h"

namespace rcspp::preprocess::ng {

struct DssrOptions {
    int32_t initial_ng_size = 4;
    int32_t max_ng_size = 12;
    int32_t dssr_iterations = 6;
    bool enable_simd = true;
};

struct DssrBoundsResult {
    std::vector<double> fwd;
    std::vector<double> bwd;
    int32_t ng_size = 1;
    bool elementary_path_found = false;
    double elementary_path_cost = std::numeric_limits<double>::infinity();
    std::vector<int32_t> elementary_path;
};

namespace detail {

constexpr double kTol = 1e-12;

struct NodeLabels {
    std::vector<double> cost;
    std::vector<double> demand;
    std::vector<int32_t> prev_node;
    std::vector<int32_t> prev_slot;
    std::vector<uint8_t> active;
    std::vector<uint64_t> visited;  // slot-major: [slot * words + w]
    int32_t active_count = 0;
};

struct LabelingRun {
    std::vector<NodeLabels> nodes;
    std::vector<double> best_cost;
    std::vector<int32_t> best_slot;
    int32_t root = -1;
    int32_t words = 0;
};

inline int32_t bitmap_words(int32_t n) { return (n + 63) / 64; }

inline void set_bit(uint64_t* bm, int32_t v) {
    bm[v >> 6] |= (uint64_t{1} << (v & 63));
}

inline bool test_bit(const uint64_t* bm, int32_t v) {
    return (bm[v >> 6] & (uint64_t{1} << (v & 63))) != 0;
}

inline bool is_subset_bits(const uint64_t* a, const uint64_t* b, int32_t words) {
    for (int32_t w = 0; w < words; ++w) {
        if ((a[w] & b[w]) != a[w]) return false;
    }
    return true;
}

inline uint64_t* visited_ptr(NodeLabels& store, int32_t slot, int32_t words) {
    return store.visited.data() + static_cast<size_t>(slot) * words;
}

inline const uint64_t* visited_ptr(const NodeLabels& store, int32_t slot, int32_t words) {
    return store.visited.data() + static_cast<size_t>(slot) * words;
}

inline void collect_candidates_leq(const NodeLabels& store,
                                   double cost_limit,
                                   double demand_limit,
                                   bool enable_simd,
                                   std::vector<int32_t>& out) {
    out.clear();
    const size_t n = store.cost.size();
#ifdef __AVX2__
    if (enable_simd) {
        const __m256d c_lim = _mm256_set1_pd(cost_limit + kTol);
        const __m256d d_lim = _mm256_set1_pd(demand_limit + kTol);
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            const __m256d c = _mm256_loadu_pd(&store.cost[i]);
            const __m256d d = _mm256_loadu_pd(&store.demand[i]);
            const __m256d c_ok = _mm256_cmp_pd(c, c_lim, _CMP_LE_OQ);
            const __m256d d_ok = _mm256_cmp_pd(d, d_lim, _CMP_LE_OQ);
            int mask = _mm256_movemask_pd(_mm256_and_pd(c_ok, d_ok));
            while (mask) {
                const int b = __builtin_ctz(mask);
                const size_t idx = i + static_cast<size_t>(b);
                if (store.active[idx]) out.push_back(static_cast<int32_t>(idx));
                mask &= (mask - 1);
            }
        }
        for (; i < n; ++i) {
            if (!store.active[i]) continue;
            if (store.cost[i] <= cost_limit + kTol && store.demand[i] <= demand_limit + kTol) {
                out.push_back(static_cast<int32_t>(i));
            }
        }
        return;
    }
#else
    (void)enable_simd;
#endif
    for (size_t i = 0; i < n; ++i) {
        if (!store.active[i]) continue;
        if (store.cost[i] <= cost_limit + kTol && store.demand[i] <= demand_limit + kTol) {
            out.push_back(static_cast<int32_t>(i));
        }
    }
}

inline void collect_candidates_geq(const NodeLabels& store,
                                   double cost_floor,
                                   double demand_floor,
                                   bool enable_simd,
                                   std::vector<int32_t>& out) {
    out.clear();
    const size_t n = store.cost.size();
#ifdef __AVX2__
    if (enable_simd) {
        const __m256d c_floor = _mm256_set1_pd(cost_floor - kTol);
        const __m256d d_floor = _mm256_set1_pd(demand_floor - kTol);
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            const __m256d c = _mm256_loadu_pd(&store.cost[i]);
            const __m256d d = _mm256_loadu_pd(&store.demand[i]);
            const __m256d c_ok = _mm256_cmp_pd(c, c_floor, _CMP_GE_OQ);
            const __m256d d_ok = _mm256_cmp_pd(d, d_floor, _CMP_GE_OQ);
            int mask = _mm256_movemask_pd(_mm256_and_pd(c_ok, d_ok));
            while (mask) {
                const int b = __builtin_ctz(mask);
                const size_t idx = i + static_cast<size_t>(b);
                if (store.active[idx]) out.push_back(static_cast<int32_t>(idx));
                mask &= (mask - 1);
            }
        }
        for (; i < n; ++i) {
            if (!store.active[i]) continue;
            if (store.cost[i] + kTol >= cost_floor && store.demand[i] + kTol >= demand_floor) {
                out.push_back(static_cast<int32_t>(i));
            }
        }
        return;
    }
#else
    (void)enable_simd;
#endif
    for (size_t i = 0; i < n; ++i) {
        if (!store.active[i]) continue;
        if (store.cost[i] + kTol >= cost_floor && store.demand[i] + kTol >= demand_floor) {
            out.push_back(static_cast<int32_t>(i));
        }
    }
}

inline std::vector<std::vector<int32_t>> nearest_neighbors(const Problem& prob) {
    const int32_t n = prob.num_nodes();
    const auto& graph = prob.graph();
    std::vector<std::vector<int32_t>> nearest(n);
    for (int32_t i = 0; i < n; ++i) {
        std::vector<std::pair<double, int32_t>> tmp;
        tmp.reserve(graph.incident_edges(i).size());
        for (auto e : graph.incident_edges(i)) {
            const int32_t j = graph.other_endpoint(e, i);
            tmp.emplace_back(prob.edge_cost(e), j);
        }
        std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
        nearest[i].reserve(tmp.size());
        for (const auto& [_, j] : tmp) nearest[i].push_back(j);
    }
    return nearest;
}

inline void initialize_ng_masks(std::vector<uint64_t>& masks,
                                int32_t n,
                                int32_t words,
                                int32_t ng_size,
                                const std::vector<std::vector<int32_t>>& nearest) {
    masks.assign(static_cast<size_t>(n) * words, 0);
    const int32_t k = std::max<int32_t>(1, ng_size);
    for (int32_t i = 0; i < n; ++i) {
        uint64_t* row = masks.data() + static_cast<size_t>(i) * words;
        set_bit(row, i);
        const int32_t add = std::max<int32_t>(0, k - 1);
        for (int32_t t = 0; t < add && t < static_cast<int32_t>(nearest[i].size()); ++t) {
            set_bit(row, nearest[i][t]);
        }
    }
}

inline int32_t row_bitcount(const std::vector<uint64_t>& masks, int32_t row, int32_t words) {
    const uint64_t* p = masks.data() + static_cast<size_t>(row) * words;
    int32_t cnt = 0;
    for (int32_t w = 0; w < words; ++w) cnt += static_cast<int32_t>(__builtin_popcountll(p[w]));
    return cnt;
}

inline int32_t max_ng_cardinality(const std::vector<uint64_t>& masks, int32_t n, int32_t words) {
    int32_t best = 0;
    for (int32_t i = 0; i < n; ++i) best = std::max(best, row_bitcount(masks, i, words));
    return best;
}

inline LabelingRun run_labeling(const Problem& prob,
                                int32_t root,
                                std::span<const double> edge_costs,
                                std::span<const double> profits,
                                const std::vector<uint64_t>& ng_masks,
                                bool enable_simd) {
    const int32_t n = prob.num_nodes();
    const double Q = prob.capacity();
    const int32_t words = bitmap_words(n);
    const auto& graph = prob.graph();

    LabelingRun run;
    run.nodes.resize(n);
    run.best_cost.assign(n, std::numeric_limits<double>::infinity());
    run.best_slot.assign(n, -1);
    run.root = root;
    run.words = words;

    auto& root_store = run.nodes[root];
    root_store.cost.push_back(-profits[root]);
    root_store.demand.push_back(prob.demand(root));
    root_store.prev_node.push_back(-1);
    root_store.prev_slot.push_back(-1);
    root_store.active.push_back(1);
    root_store.visited.assign(words, 0);
    set_bit(root_store.visited.data(), root);
    root_store.active_count = 1;
    run.best_cost[root] = root_store.cost[0];
    run.best_slot[root] = 0;

    struct QueueEntry {
        int32_t node;
        int32_t slot;
    };
    std::vector<QueueEntry> queue;
    queue.push_back({root, 0});

    std::vector<int32_t> candidate_idx;
    std::vector<uint64_t> temp_vis(words, 0);
    size_t head = 0;
    while (head < queue.size()) {
        const auto [u, slot_u] = queue[head++];
        auto& src_store = run.nodes[u];
        if (slot_u < 0 || slot_u >= static_cast<int32_t>(src_store.cost.size())) continue;
        if (!src_store.active[slot_u]) continue;

        const double base_cost = src_store.cost[slot_u];
        const double base_dem = src_store.demand[slot_u];
        const int32_t prev_u = src_store.prev_node[slot_u];
        const uint64_t* vis_u = visited_ptr(src_store, slot_u, words);

        for (auto e : graph.incident_edges(u)) {
            const int32_t v = graph.other_endpoint(e, u);
            if (v == prev_u) continue;
            const double new_demand = base_dem + prob.demand(v);
            if (new_demand > Q) continue;
            if (test_bit(vis_u, v)) continue;

            const double new_cost = base_cost + edge_costs[e] - profits[v];
            const uint64_t* ng_row = ng_masks.data() + static_cast<size_t>(v) * words;
            for (int32_t w = 0; w < words; ++w) temp_vis[w] = vis_u[w] & ng_row[w];
            set_bit(temp_vis.data(), v);

            auto& dst_store = run.nodes[v];
            bool dominated = false;
            collect_candidates_leq(dst_store, new_cost, new_demand, enable_simd, candidate_idx);
            for (int32_t idx : candidate_idx) {
                const uint64_t* ex_vis = visited_ptr(dst_store, idx, words);
                if (is_subset_bits(ex_vis, temp_vis.data(), words)) {
                    dominated = true;
                    break;
                }
            }
            if (dominated) continue;

            collect_candidates_geq(dst_store, new_cost, new_demand, enable_simd, candidate_idx);
            for (int32_t idx : candidate_idx) {
                const uint64_t* ex_vis = visited_ptr(dst_store, idx, words);
                if (is_subset_bits(temp_vis.data(), ex_vis, words)) {
                    if (dst_store.active[idx]) {
                        dst_store.active[idx] = 0;
                        dst_store.active_count--;
                    }
                }
            }

            const int32_t new_slot = static_cast<int32_t>(dst_store.cost.size());
            dst_store.cost.push_back(new_cost);
            dst_store.demand.push_back(new_demand);
            dst_store.prev_node.push_back(u);
            dst_store.prev_slot.push_back(slot_u);
            dst_store.active.push_back(1);
            dst_store.visited.resize(dst_store.visited.size() + static_cast<size_t>(words));
            std::copy(temp_vis.begin(), temp_vis.end(), dst_store.visited.begin() +
                                                     static_cast<ptrdiff_t>(new_slot) * words);
            dst_store.active_count++;

            if (new_cost < run.best_cost[v] - kTol) {
                run.best_cost[v] = new_cost;
                run.best_slot[v] = new_slot;
            }
            queue.push_back({v, new_slot});
        }
    }

    return run;
}

inline std::vector<int32_t> reconstruct_best_path(const LabelingRun& run, int32_t node) {
    std::vector<int32_t> rev;
    if (node < 0 || node >= static_cast<int32_t>(run.best_slot.size())) return rev;
    int32_t slot = run.best_slot[node];
    if (slot < 0) return rev;
    int32_t cur = node;
    int32_t guard = 0;
    const int32_t guard_limit = std::max<int32_t>(1, static_cast<int32_t>(run.best_slot.size()) * 16);
    while (cur >= 0 && slot >= 0 && guard++ < guard_limit) {
        rev.push_back(cur);
        const auto& store = run.nodes[cur];
        const int32_t next_node = store.prev_node[slot];
        const int32_t next_slot = store.prev_slot[slot];
        cur = next_node;
        slot = next_slot;
    }
    std::reverse(rev.begin(), rev.end());
    return rev;
}

inline bool is_elementary_path(const std::vector<int32_t>& path) {
    std::unordered_map<int32_t, uint8_t> seen;
    seen.reserve(path.size() * 2 + 1);
    for (int32_t v : path) {
        if (seen.contains(v)) return false;
        seen.emplace(v, 1);
    }
    return true;
}

inline void collect_cycles(const LabelingRun& run, std::vector<std::vector<int32_t>>& cycles) {
    const int32_t n = static_cast<int32_t>(run.best_slot.size());
    cycles.clear();
    for (int32_t v = 0; v < n; ++v) {
        if (run.best_slot[v] < 0) continue;
        auto path = reconstruct_best_path(run, v);
        if (path.size() < 3) continue;
        std::unordered_map<int32_t, int32_t> seen;
        seen.reserve(path.size() * 2 + 1);
        for (int32_t i = 0; i < static_cast<int32_t>(path.size()); ++i) {
            auto it = seen.find(path[i]);
            if (it != seen.end()) {
                std::vector<int32_t> cyc;
                cyc.reserve(static_cast<size_t>(i - it->second + 1));
                for (int32_t j = it->second; j <= i; ++j) cyc.push_back(path[j]);
                if (cyc.size() >= 2) cycles.push_back(std::move(cyc));
                break;
            }
            seen.emplace(path[i], i);
        }
    }
}

inline bool grow_ng_from_cycles(std::vector<uint64_t>& ng_masks,
                                int32_t n,
                                int32_t words,
                                const std::vector<std::vector<int32_t>>& nearest,
                                const std::vector<std::vector<int32_t>>& cycles,
                                int32_t max_ng_size) {
    bool changed = false;
    for (const auto& cycle : cycles) {
        if (cycle.size() < 2) continue;
        for (int32_t u : cycle) {
            if (u < 0 || u >= n) continue;
            uint64_t* row = ng_masks.data() + static_cast<size_t>(u) * words;
            int32_t cnt = row_bitcount(ng_masks, u, words);
            for (int32_t w : cycle) {
                if (cnt >= max_ng_size) break;
                if (!test_bit(row, w)) {
                    set_bit(row, w);
                    cnt++;
                    changed = true;
                }
            }
            for (int32_t nb : nearest[u]) {
                if (cnt >= max_ng_size) break;
                if (!test_bit(row, nb)) {
                    set_bit(row, nb);
                    cnt++;
                    changed = true;
                }
            }
        }
    }
    return changed;
}

}  // namespace detail

inline DssrBoundsResult compute_bounds(const Problem& prob,
                                       int32_t source,
                                       int32_t target,
                                       std::span<const double> edge_costs,
                                       std::span<const double> profits,
                                       const DssrOptions& opts) {
    const int32_t n = prob.num_nodes();
    const int32_t words = detail::bitmap_words(n);
    const auto nearest = detail::nearest_neighbors(prob);

    const int32_t init_ng = std::max<int32_t>(1, std::min(opts.initial_ng_size, opts.max_ng_size));
    const int32_t max_ng = std::max<int32_t>(init_ng, opts.max_ng_size);
    const int32_t max_iter = std::max<int32_t>(1, opts.dssr_iterations);

    std::vector<uint64_t> ng_masks;
    detail::initialize_ng_masks(ng_masks, n, words, init_ng, nearest);

    DssrBoundsResult result;
    result.ng_size = detail::max_ng_cardinality(ng_masks, n, words);

    std::vector<std::vector<int32_t>> cycles;

    for (int32_t it = 0; it < max_iter; ++it) {
        auto fwd_run = detail::run_labeling(
            prob, source, edge_costs, profits, ng_masks, opts.enable_simd);
        result.fwd = fwd_run.best_cost;

        if (source == target) {
            result.bwd = result.fwd;
        } else {
            auto bwd_run = detail::run_labeling(
                prob, target, edge_costs, profits, ng_masks, opts.enable_simd);
            result.bwd = bwd_run.best_cost;
            if (fwd_run.best_slot[target] >= 0) {
                auto path = detail::reconstruct_best_path(fwd_run, target);
                if (!path.empty() && detail::is_elementary_path(path)) {
                    const double path_cost = fwd_run.best_cost[target];
                    if (!result.elementary_path_found
                        || path_cost < result.elementary_path_cost - detail::kTol) {
                        result.elementary_path_found = true;
                        result.elementary_path_cost = path_cost;
                        result.elementary_path = std::move(path);
                    }
                }
            }
            if (it + 1 >= max_iter) break;
            detail::collect_cycles(fwd_run, cycles);
            std::vector<std::vector<int32_t>> cycles_bwd;
            detail::collect_cycles(bwd_run, cycles_bwd);
            cycles.insert(cycles.end(),
                          std::make_move_iterator(cycles_bwd.begin()),
                          std::make_move_iterator(cycles_bwd.end()));
        }

        if (source == target && it + 1 >= max_iter) break;
        if (source == target) {
            detail::collect_cycles(fwd_run, cycles);
        }
        if (it + 1 >= max_iter) break;
        if (cycles.empty()) break;

        bool grown = detail::grow_ng_from_cycles(
            ng_masks, n, words, nearest, cycles, max_ng);
        result.ng_size = detail::max_ng_cardinality(ng_masks, n, words);
        if (!grown || result.ng_size >= max_ng) break;
    }

    if (result.fwd.empty()) {
        auto fwd_run = detail::run_labeling(
            prob, source, edge_costs, profits, ng_masks, opts.enable_simd);
        result.fwd = std::move(fwd_run.best_cost);
    }
    if (result.bwd.empty()) {
        if (source == target) {
            result.bwd = result.fwd;
        } else {
            auto bwd_run = detail::run_labeling(
                prob, target, edge_costs, profits, ng_masks, opts.enable_simd);
            result.bwd = std::move(bwd_run.best_cost);
        }
    }
    return result;
}

inline DssrBoundsResult compute_bounds(const Problem& prob,
                                       int32_t source,
                                       int32_t target,
                                       const DssrOptions& opts) {
    return compute_bounds(prob, source, target, prob.edge_costs(), prob.profits(), opts);
}

}  // namespace rcspp::preprocess::ng
