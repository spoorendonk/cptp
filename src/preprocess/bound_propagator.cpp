#include "preprocess/bound_propagator.h"

#include <algorithm>
#include <cstdint>

namespace cptp::preprocess {

BoundPropagator::BoundPropagator(const Problem& prob,
                                 std::vector<double> fwd_bounds,
                                 std::vector<double> bwd_bounds,
                                 double correction)
    : prob_(prob),
      fwd_bounds_(std::move(fwd_bounds)),
      bwd_bounds_(std::move(bwd_bounds)),
      correction_(correction) {
  const auto& graph = prob_.graph();
  const int32_t m = prob_.num_edges();
  const int32_t n = prob_.num_nodes();

  // Pre-build adjacency
  adjacency_.resize(n);
  for (auto e : graph.edges()) {
    int32_t u = graph.edge_source(e);
    int32_t v = graph.edge_target(e);
    adjacency_[u].push_back({e, v});
    adjacency_[v].push_back({e, u});
  }

  // Pre-build edge costs
  edge_costs_.resize(m);
  for (auto e : graph.edges()) {
    edge_costs_[e] = prob_.edge_cost(e);
  }

  // Pre-build profits
  profits_.resize(n);
  for (int32_t i = 0; i < n; ++i) {
    profits_[i] = prob_.profit(i);
  }
}

void BoundPropagator::set_all_pairs_bounds(std::vector<double> dist) {
  all_pairs_ = std::move(dist);
}

std::vector<int32_t> BoundPropagator::sweep(
    double upper_bound, std::span<const double> col_upper) const {
  const int32_t m = prob_.num_edges();
  const auto& f = fwd_bounds_;
  const auto& b = bwd_bounds_;

  std::vector<int32_t> fixings;

  for (int32_t e = 0; e < m; ++e) {
    if (col_upper[e] < 0.5) continue;  // already fixed

    int32_t u = prob_.graph().edge_source(e);
    int32_t v = prob_.graph().edge_target(e);

    double lb1 = f[u] + edge_costs_[e] + b[v] + correction_;
    double lb2 = f[v] + edge_costs_[e] + b[u] + correction_;
    double lb = std::min(lb1, lb2);

    if (lb > upper_bound + 1e-6) {
      fixings.push_back(e);
    }
  }

  return fixings;
}

std::vector<int32_t> BoundPropagator::sweep_nodes(
    std::span<const double> col_upper, int32_t y_offset) const {
  const int32_t n = prob_.num_nodes();
  std::vector<int32_t> fixings;

  for (int32_t i = 0; i < n; ++i) {
    if (col_upper[y_offset + i] < 0.5) continue;
    bool all_fixed = true;
    for (const auto& [e, nb] : adjacency_[i]) {
      if (col_upper[e] > 0.5) {
        all_fixed = false;
        break;
      }
    }
    if (all_fixed) {
      fixings.push_back(y_offset + i);
    }
  }

  return fixings;
}

std::vector<int32_t> BoundPropagator::propagate_fixed_edge(
    int32_t edge, double upper_bound, std::span<const double> col_upper) const {
  int32_t a = prob_.graph().edge_source(edge);
  int32_t i = prob_.graph().edge_target(edge);

  std::vector<int32_t> fixings;

  if (!all_pairs_.empty()) {
    // All-pairs Trigger B: scan ALL unfixed edges
    const int32_t m = prob_.num_edges();
    const int32_t n = prob_.num_nodes();
    const int32_t depot = prob_.depot();
    const double d_depot_a = all_pairs_[depot * n + a];
    const double d_depot_i = all_pairs_[depot * n + i];
    const double d_i_depot = all_pairs_[i * n + depot];
    const double d_a_depot = all_pairs_[a * n + depot];
    const double cost_ai = edge_costs_[edge];

    for (int32_t ej = 0; ej < m; ++ej) {
      if (ej == edge) continue;
      if (col_upper[ej] < 0.5) continue;

      int32_t u = prob_.graph().edge_source(ej);
      int32_t v = prob_.graph().edge_target(ej);

      double lb1 = d_depot_a + cost_ai + all_pairs_[i * n + u] +
                   edge_costs_[ej] + all_pairs_[v * n + depot] + correction_;
      double lb2 = all_pairs_[depot * n + u] + edge_costs_[ej] +
                   all_pairs_[v * n + a] + cost_ai + d_i_depot + correction_;
      double lb3 = d_depot_a + cost_ai + all_pairs_[i * n + v] +
                   edge_costs_[ej] + all_pairs_[u * n + depot] + correction_;
      double lb4 = all_pairs_[depot * n + v] + edge_costs_[ej] +
                   all_pairs_[u * n + a] + cost_ai + d_i_depot + correction_;
      double lb5 = d_depot_i + cost_ai + all_pairs_[a * n + u] +
                   edge_costs_[ej] + all_pairs_[v * n + depot] + correction_;
      double lb6 = all_pairs_[depot * n + u] + edge_costs_[ej] +
                   all_pairs_[v * n + i] + cost_ai + d_a_depot + correction_;
      double lb7 = d_depot_i + cost_ai + all_pairs_[a * n + v] +
                   edge_costs_[ej] + all_pairs_[u * n + depot] + correction_;
      double lb8 = all_pairs_[depot * n + v] + edge_costs_[ej] +
                   all_pairs_[u * n + i] + cost_ai + d_a_depot + correction_;

      double lb = std::min({lb1, lb2, lb3, lb4, lb5, lb6, lb7, lb8});

      if (lb > upper_bound + 1e-6) {
        fixings.push_back(ej);
      }
    }
  } else {
    // Fallback: neighbor-only scan (original Trigger B)
    const auto& f = fwd_bounds_;
    const auto& b = bwd_bounds_;
    double cost_a_to_i = f[a] + edge_costs_[edge] - profits_[i];

    for (const auto& [ej, j] : adjacency_[i]) {
      if (ej == edge) continue;
      if (col_upper[ej] < 0.5) continue;

      double lb = cost_a_to_i + edge_costs_[ej] + b[j] + correction_;
      if (lb > upper_bound + 1e-6) {
        fixings.push_back(ej);
      }
    }

    double cost_via_i_return = edge_costs_[edge] + b[i] - profits_[a];

    for (const auto& [ek, k] : adjacency_[a]) {
      if (ek == edge) continue;
      if (col_upper[ek] < 0.5) continue;

      double lb = f[k] + edge_costs_[ek] + cost_via_i_return + correction_;
      if (lb > upper_bound + 1e-6) {
        fixings.push_back(ek);
      }
    }
  }

  return fixings;
}

}  // namespace cptp::preprocess
