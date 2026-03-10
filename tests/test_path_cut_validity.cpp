#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "core/digraph.h"
#include "core/gomory_hu.h"
#include "core/problem.h"
#include "preprocess/edge_elimination.h"
#include "sep/comb_separator.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/rglm_separator.h"
#include "sep/sec_separator.h"
#include "sep/separation_context.h"
#include "sep/spi_separator.h"

namespace {

struct SupportData {
  cptp::digraph graph;
  std::vector<double> capacity;
  std::unique_ptr<cptp::gomory_hu_tree> tree;
};

struct IntegerPathSolution {
  std::vector<double> x_values;
  std::vector<double> y_values;
};

struct EquivalentInstances {
  cptp::Problem tour;
  cptp::Problem path;
};

struct EquivalentSample {
  std::vector<double> x_tour;
  std::vector<double> y_tour;
  std::vector<double> x_path;
  std::vector<double> y_path;
};

SupportData build_support(const cptp::Problem& prob,
                          std::span<const double> x_values, double tol = 1e-6) {
  const auto& g = prob.graph();
  cptp::digraph_builder builder(prob.num_nodes());
  for (auto e : g.edges()) {
    const double xval = x_values[e];
    if (xval <= tol) continue;
    const int32_t u = g.edge_source(e);
    const int32_t v = g.edge_target(e);
    builder.add_arc(u, v, xval);
    builder.add_arc(v, u, xval);
  }
  auto [support_graph, cap] = builder.build();
  auto tree =
      std::make_unique<cptp::gomory_hu_tree>(support_graph, cap, prob.source());
  return {std::move(support_graph), std::move(cap), std::move(tree)};
}

double evaluate_cut_lhs(const cptp::sep::Cut& cut,
                        std::span<const double> x_values,
                        std::span<const double> y_values) {
  const int32_t m = static_cast<int32_t>(x_values.size());
  double lhs = 0.0;
  for (size_t k = 0; k < cut.indices.size(); ++k) {
    const int32_t idx = cut.indices[k];
    const double coeff = cut.values[k];
    lhs += (idx < m) ? coeff * x_values[idx] : coeff * y_values[idx - m];
  }
  return lhs;
}

int64_t quantize_coeff(double v) {
  return static_cast<int64_t>(std::llround(v * 1e9));
}

using CutSignature =
    std::pair<int64_t, std::vector<std::pair<int32_t, int64_t>>>;

std::vector<CutSignature> normalize_cut_signatures(
    const std::vector<cptp::sep::Cut>& cuts) {
  std::vector<CutSignature> out;
  out.reserve(cuts.size());
  for (const auto& cut : cuts) {
    std::vector<std::pair<int32_t, int64_t>> terms;
    terms.reserve(cut.indices.size());
    for (size_t k = 0; k < cut.indices.size(); ++k) {
      terms.push_back({cut.indices[k], quantize_coeff(cut.values[k])});
    }
    std::sort(terms.begin(), terms.end());

    // Merge duplicates if any separator emits repeated variable indices.
    std::vector<std::pair<int32_t, int64_t>> merged;
    for (const auto& t : terms) {
      if (!merged.empty() && merged.back().first == t.first) {
        merged.back().second += t.second;
      } else {
        merged.push_back(t);
      }
    }
    out.push_back({quantize_coeff(cut.rhs), std::move(merged)});
  }
  std::sort(out.begin(), out.end());
  return out;
}

int32_t edge_index(const cptp::Problem& prob, int32_t a, int32_t b) {
  if (a > b) std::swap(a, b);
  const auto& g = prob.graph();
  for (int32_t e : g.incident_edges(a)) {
    if (g.other_endpoint(e, a) == b) return e;
  }
  return -1;
}

cptp::Problem make_validation_path_problem() {
  cptp::Problem prob;
  std::vector<cptp::Edge> edges;
  std::vector<double> costs;

  // Complete graph on 6 nodes: source=0, target=5.
  // Internal nodes 1..4 are customers.
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      edges.push_back({i, j});
      // Keep chain-ish edges cheap, shortcuts expensive.
      if (j == i + 1) {
        costs.push_back(4.0);
      } else {
        costs.push_back(15.0);
      }
    }
  }

  std::vector<double> profits = {0, 2, 2, 2, 2, 0};
  std::vector<double> demands = {0, 2, 2, 2, 2, 0};
  constexpr double capacity = 6.0;
  prob.build(6, edges, costs, profits, demands, capacity, 0, 5);
  return prob;
}

std::vector<IntegerPathSolution> enumerate_feasible_paths(
    const cptp::Problem& prob) {
  const auto& g = prob.graph();
  const int32_t n = prob.num_nodes();
  const int32_t m = prob.num_edges();
  const int32_t source = prob.source();
  const int32_t target = prob.target();
  const double cap = prob.capacity();

  std::vector<IntegerPathSolution> out;
  std::vector<bool> visited(n, false);
  std::vector<int32_t> path_nodes;
  std::vector<int32_t> path_edges;
  visited[source] = true;
  path_nodes.push_back(source);

  auto push_solution = [&]() {
    double demand_sum = 0.0;
    for (int32_t node : path_nodes) demand_sum += prob.demand(node);
    if (demand_sum > cap + 1e-9) return;

    IntegerPathSolution sol;
    sol.x_values.assign(m, 0.0);
    sol.y_values.assign(n, 0.0);
    for (int32_t e : path_edges) sol.x_values[e] = 1.0;
    for (int32_t node : path_nodes) sol.y_values[node] = 1.0;
    out.push_back(std::move(sol));
  };

  std::function<void(int32_t)> dfs = [&](int32_t u) {
    if (u == target) {
      push_solution();
      return;
    }

    for (int32_t e : g.incident_edges(u)) {
      const int32_t v = g.other_endpoint(e, u);
      if (visited[v]) continue;
      visited[v] = true;
      path_nodes.push_back(v);
      path_edges.push_back(e);
      dfs(v);
      path_edges.pop_back();
      path_nodes.pop_back();
      visited[v] = false;
    }
  };

  dfs(source);
  return out;
}

EquivalentInstances make_equivalent_instances() {
  // Tour: depot 0, customers 1..4.
  cptp::Problem tour;
  std::vector<cptp::Edge> edges_tour;
  std::vector<double> costs_tour;
  for (int i = 0; i <= 4; ++i) {
    for (int j = i + 1; j <= 4; ++j) {
      edges_tour.push_back({i, j});
      if (i == 0 || j == 0) {
        const int c = (i == 0) ? j : i;
        costs_tour.push_back(2.0 + static_cast<double>(c));
      } else {
        costs_tour.push_back(3.0 + std::abs(i - j));
      }
    }
  }
  std::vector<double> profits_tour = {0, 2, 2, 2, 2};
  std::vector<double> demands_tour = {0, 2, 2, 2, 2};
  tour.build(5, edges_tour, costs_tour, profits_tour, demands_tour, 7.0, 0, 0);

  // Path equivalent via depot split: source 0, customers 1..4, target 5.
  cptp::Problem path;
  std::vector<cptp::Edge> edges_path;
  std::vector<double> costs_path;

  for (int i = 1; i <= 4; ++i) {
    for (int j = i + 1; j <= 4; ++j) {
      edges_path.push_back({i, j});
      costs_path.push_back(3.0 + std::abs(i - j));
    }
  }
  for (int i = 1; i <= 4; ++i) {
    edges_path.push_back({0, i});  // source copy of depot-customer
    costs_path.push_back(2.0 + static_cast<double>(i));
    edges_path.push_back({i, 5});  // target copy of depot-customer
    costs_path.push_back(2.0 + static_cast<double>(i));
  }
  edges_path.push_back(
      {0, 5});  // keep terminals tightly connected in support graph
  costs_path.push_back(0.5);

  std::vector<double> profits_path = {0, 2, 2, 2, 2, 0};
  std::vector<double> demands_path = {0, 2, 2, 2, 2, 0};
  path.build(6, edges_path, costs_path, profits_path, demands_path, 7.0, 0, 5);
  return {std::move(tour), std::move(path)};
}

EquivalentInstances make_comb_equivalent_instances() {
  // Tour with 6 customers (1..6) and depot 0.
  cptp::Problem tour;
  std::vector<cptp::Edge> edges_tour;
  std::vector<double> costs_tour;
  for (int i = 0; i <= 6; ++i) {
    for (int j = i + 1; j <= 6; ++j) {
      edges_tour.push_back({i, j});
      costs_tour.push_back(1.0 + std::abs(i - j));
    }
  }
  std::vector<double> profits_tour = {0, 1, 1, 1, 1, 1, 1};
  std::vector<double> demands_tour = {0, 1, 1, 1, 1, 1, 1};
  tour.build(7, edges_tour, costs_tour, profits_tour, demands_tour, 100.0, 0,
             0);

  // Depot-split path equivalent: source 0, target 7, customers 1..6.
  cptp::Problem path;
  std::vector<cptp::Edge> edges_path;
  std::vector<double> costs_path;
  for (int i = 1; i <= 6; ++i) {
    for (int j = i + 1; j <= 6; ++j) {
      edges_path.push_back({i, j});
      costs_path.push_back(1.0 + std::abs(i - j));
    }
  }
  for (int i = 1; i <= 6; ++i) {
    edges_path.push_back({0, i});
    costs_path.push_back(1.0 + i);
    edges_path.push_back({i, 7});
    costs_path.push_back(1.0 + i);
  }
  edges_path.push_back({0, 7});
  costs_path.push_back(0.5);

  std::vector<double> profits_path = {0, 1, 1, 1, 1, 1, 1, 0};
  std::vector<double> demands_path = {0, 1, 1, 1, 1, 1, 1, 0};
  path.build(8, edges_path, costs_path, profits_path, demands_path, 100.0, 0,
             7);
  return {std::move(tour), std::move(path)};
}

EquivalentSample make_comb_equivalent_sample(const EquivalentInstances& inst) {
  EquivalentSample sample;
  sample.x_tour.assign(inst.tour.num_edges(), 0.02);
  sample.y_tour.assign(inst.tour.num_nodes(), 0.2);
  sample.x_path.assign(inst.path.num_edges(), 0.02);
  sample.y_path.assign(inst.path.num_nodes(), 0.2);

  sample.y_tour[0] = 1.0;
  sample.y_path[0] = 1.0;
  sample.y_path[inst.path.target()] = 1.0;

  auto set_tour = [&](int a, int b, double v) {
    const int32_t e = edge_index(inst.tour, a, b);
    REQUIRE(e >= 0);
    sample.x_tour[e] = v;
  };
  auto set_path = [&](int a, int b, double v) {
    const int32_t e = edge_index(inst.path, a, b);
    REQUIRE(e >= 0);
    sample.x_path[e] = v;
  };

  // Handle H = {terminal, 1,2,3}. Keep depot-split mapping exact:
  // x_tour(0,i) = x_path(0,i) + x_path(i,7).
  for (int u : {1, 2, 3}) {
    set_tour(0, u, 0.9);
    set_path(0, u, 0.45);
    set_path(u, 7, 0.45);
  }
  // Keep outside customers out of the handle while preserving split mapping.
  for (int u : {4, 5, 6}) {
    set_tour(0, u, 0.05);
    set_path(0, u, 0.025);
    set_path(u, 7, 0.025);
  }
  set_path(0, 7, 0.0);

  // Dense inside structure boosts x(E(H)).
  for (auto [a, b] : std::vector<std::pair<int, int>>{{1, 2}, {1, 3}, {2, 3}}) {
    set_tour(a, b, 0.85);
    set_path(a, b, 0.85);
  }

  // Three disjoint teeth to outside customers.
  // Keep them below Comb BFS thresholds (0.3/0.5) so outside nodes stay outside
  // H.
  set_tour(1, 4, 0.25);
  set_path(1, 4, 0.25);
  set_tour(2, 5, 0.25);
  set_path(2, 5, 0.25);
  set_tour(3, 6, 0.25);
  set_path(3, 6, 0.25);

  return sample;
}

std::vector<EquivalentSample> build_equivalent_samples(
    const EquivalentInstances& inst, int num_samples) {
  const auto& tour = inst.tour;
  const auto& path = inst.path;
  const int32_t m_tour = tour.num_edges();
  const int32_t m_path = path.num_edges();

  std::vector<int32_t> idx_t_depot(5, -1);
  std::vector<int32_t> idx_p_src(5, -1), idx_p_tgt(5, -1);
  int32_t idx_p_st = edge_index(path, 0, 5);
  int32_t idx_t_cc[5][5];
  int32_t idx_p_cc[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      idx_t_cc[i][j] = -1;
      idx_p_cc[i][j] = -1;
    }
  }

  for (int c = 1; c <= 4; ++c) {
    idx_t_depot[c] = edge_index(tour, 0, c);
    idx_p_src[c] = edge_index(path, 0, c);
    idx_p_tgt[c] = edge_index(path, c, 5);
  }
  for (int i = 1; i <= 4; ++i) {
    for (int j = i + 1; j <= 4; ++j) {
      idx_t_cc[i][j] = edge_index(tour, i, j);
      idx_p_cc[i][j] = edge_index(path, i, j);
    }
  }

  std::mt19937 rng(424242);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::vector<EquivalentSample> out;
  out.reserve(static_cast<size_t>(num_samples));

  for (int s = 0; s < num_samples; ++s) {
    EquivalentSample sample;
    sample.x_tour.assign(m_tour, 0.0);
    sample.y_tour.assign(5, 0.0);
    sample.x_path.assign(m_path, 0.0);
    sample.y_path.assign(6, 0.0);

    sample.y_tour[0] = 1.0;
    sample.y_path[0] = 1.0;
    sample.y_path[5] = 1.0;
    for (int c = 1; c <= 4; ++c) {
      const double y = 0.2 + 0.75 * u01(rng);
      sample.y_tour[c] = y;
      sample.y_path[c] = y;
    }

    // Depot/customer aggregate capacities shared across both models.
    for (int c = 1; c <= 4; ++c) {
      const double a = 0.2 + 0.8 * u01(rng);
      sample.x_tour[idx_t_depot[c]] = a;
      sample.x_path[idx_p_src[c]] = 0.5 * a;
      sample.x_path[idx_p_tgt[c]] = 0.5 * a;
    }
    for (int i = 1; i <= 4; ++i) {
      for (int j = i + 1; j <= 4; ++j) {
        const double b =
            (u01(rng) < 0.6) ? (0.55 + 0.45 * u01(rng)) : (0.02 * u01(rng));
        sample.x_tour[idx_t_cc[i][j]] = b;
        sample.x_path[idx_p_cc[i][j]] = b;
      }
    }
    if (idx_p_st >= 0) sample.x_path[idx_p_st] = 1.0;

    out.push_back(std::move(sample));
  }

  return out;
}

struct SetStats {
  double dS = 0.0;
  double served_demand = 0.0;
  double cut_flow = 0.0;
  double star_sum = 0.0;
};

SetStats compute_set_stats(const cptp::Problem& prob,
                           std::span<const double> x_values,
                           std::span<const double> y_values,
                           std::span<const uint8_t> in_S) {
  SetStats s;
  const auto& g = prob.graph();
  const int32_t n = prob.num_nodes();
  for (int32_t i = 0; i < n; ++i) {
    if (!in_S[i]) continue;
    const double d = prob.demand(i);
    s.dS += d;
    s.served_demand += d * y_values[i];
  }

  for (auto e : g.edges()) {
    const int32_t u = g.edge_source(e);
    const int32_t v = g.edge_target(e);
    const bool u_in = in_S[u];
    const bool v_in = in_S[v];
    if (u_in == v_in) continue;
    const int32_t outside = u_in ? v : u;
    const double x = x_values[e];
    s.cut_flow += x;
    s.star_sum += prob.demand(outside) * x;
  }
  return s;
}

std::vector<std::pair<std::vector<double>, std::vector<double>>>
build_fractional_samples(const cptp::Problem& prob) {
  const int32_t m = prob.num_edges();
  const int32_t n = prob.num_nodes();
  const int32_t source = prob.source();
  const int32_t target = prob.target();
  std::vector<std::pair<std::vector<double>, std::vector<double>>> samples;

  // Crafted sample to trigger comb-like handle/teeth structure.
  {
    std::vector<double> x(m, 0.0);
    std::vector<double> y(n, 0.0);
    y[source] = 1.0;
    y[target] = 1.0;
    for (int32_t i = 1; i < n - 1; ++i) y[i] = 0.6;

    const auto& g = prob.graph();
    for (auto e : g.edges()) {
      const int32_t u = g.edge_source(e);
      const int32_t v = g.edge_target(e);
      if ((u == 0 && (v == 1 || v == 2 || v == 3)) ||
          (v == 0 && (u == 1 || u == 2 || u == 3))) {
        x[e] = 0.85;
      } else if ((u == 1 && v == 4) || (u == 2 && v == 4) ||
                 (u == 3 && v == 4) || (u == 1 && v == 5) ||
                 (u == 2 && v == 5) || (u == 3 && v == 5)) {
        x[e] = 0.75;
      } else {
        x[e] = 0.08;
      }
    }
    samples.push_back({std::move(x), std::move(y)});
  }

  std::mt19937 rng(1234567);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  for (int sample = 0; sample < 120; ++sample) {
    std::vector<double> x(m, 0.0);
    std::vector<double> y(n, 0.0);
    for (int32_t i = 0; i < n; ++i) y[i] = 0.1 + 0.8 * u01(rng);
    y[source] = 1.0;
    y[target] = 1.0;

    // Sparse/high mixture so cuts are likely to appear.
    for (int32_t e = 0; e < m; ++e) {
      const double p = u01(rng);
      x[e] = (p < 0.45) ? (0.55 + 0.45 * u01(rng)) : (0.03 * u01(rng));
    }
    samples.push_back({std::move(x), std::move(y)});
  }

  return samples;
}

std::vector<double> compute_all_pairs(const cptp::Problem& prob) {
  const int32_t n = prob.num_nodes();
  constexpr double inf = std::numeric_limits<double>::infinity();
  std::vector<double> all_pairs(static_cast<size_t>(n) * n, inf);
  for (int32_t s = 0; s < n; ++s) {
    auto row = cptp::preprocess::forward_labeling(prob, s);
    std::copy(row.begin(), row.end(),
              all_pairs.begin() + static_cast<ptrdiff_t>(s) * n);
  }
  return all_pairs;
}

double compute_comb_lhs(const cptp::Problem& prob, std::span<const double> x,
                        std::span<const double> y,
                        std::span<const uint8_t> in_handle,
                        std::span<const std::pair<int32_t, int32_t>> teeth) {
  double lhs = 0.0;
  const auto& g = prob.graph();
  const int32_t n = prob.num_nodes();

  for (auto e : g.edges()) {
    const int32_t u = g.edge_source(e);
    const int32_t v = g.edge_target(e);
    if (in_handle[u] && in_handle[v]) lhs += x[e];
  }
  for (const auto& [a, b] : teeth) {
    const int32_t e = edge_index(prob, a, b);
    REQUIRE(e >= 0);
    lhs += x[e];
  }
  for (int32_t i = 0; i < n; ++i) {
    if (in_handle[i]) lhs -= y[i];
  }
  for (const auto& [a, b] : teeth) {
    const bool a_in = in_handle[a] != 0;
    const bool b_in = in_handle[b] != 0;
    REQUIRE(a_in != b_in);
    const int32_t out = a_in ? b : a;
    lhs -= y[out];
  }
  return lhs;
}

double compute_comb_lhs_path_split(
    const cptp::Problem& prob, std::span<const double> x,
    std::span<const double> y, std::span<const uint8_t> in_handle,
    std::span<const std::pair<int32_t, int32_t>> teeth) {
  double lhs = 0.0;
  const auto& g = prob.graph();
  const int32_t n = prob.num_nodes();
  const int32_t source = prob.source();
  const int32_t target = prob.target();

  for (auto e : g.edges()) {
    const int32_t u = g.edge_source(e);
    const int32_t v = g.edge_target(e);
    if (!in_handle[u] || !in_handle[v]) continue;
    if ((u == source && v == target) || (u == target && v == source)) continue;
    lhs += x[e];
  }
  for (const auto& [a, b] : teeth) {
    const int32_t e = edge_index(prob, a, b);
    REQUIRE(e >= 0);
    lhs += x[e];
  }
  for (int32_t i = 0; i < n; ++i) {
    if (i == target || !in_handle[i]) continue;
    lhs -= y[i];
  }
  for (const auto& [a, b] : teeth) {
    const bool a_in = in_handle[a] != 0;
    const bool b_in = in_handle[b] != 0;
    REQUIRE(a_in != b_in);
    const int32_t out = a_in ? b : a;
    REQUIRE(out != target);
    lhs -= y[out];
  }
  return lhs;
}

template <typename SeparatorT>
void validate_separator_on_path(
    const cptp::Problem& prob,
    const std::vector<IntegerPathSolution>& feasible_paths,
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>&
        samples,
    SeparatorT separator, bool require_generated_cuts = true) {
  const int32_t m = prob.num_edges();
  size_t produced_cuts = 0;

  for (const auto& [x_values, y_values] : samples) {
    auto support = build_support(prob, x_values);
    cptp::sep::SeparationContext ctx{
        .problem = prob,
        .x_values = x_values,
        .y_values = y_values,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };

    const auto cuts = separator.separate(ctx);
    produced_cuts += cuts.size();
    for (const auto& cut : cuts) {
      for (const auto& path_sol : feasible_paths) {
        const double lhs =
            evaluate_cut_lhs(cut, path_sol.x_values, path_sol.y_values);
        REQUIRE(lhs <= cut.rhs + 1e-6);
      }
    }
  }

  if (require_generated_cuts) {
    REQUIRE(produced_cuts > 0);
  }
}

}  // namespace

TEST_CASE(
    "Path cut families generate only valid cuts on all feasible s-t paths",
    "[separators][path][validity]") {
  auto prob = make_validation_path_problem();
  REQUIRE_FALSE(prob.is_tour());
  auto feasible_paths = enumerate_feasible_paths(prob);
  REQUIRE_FALSE(feasible_paths.empty());
  auto samples = build_fractional_samples(prob);

  SECTION("RCI") {
    validate_separator_on_path(prob, feasible_paths, samples,
                               cptp::sep::RCISeparator{});
  }
  SECTION("Multistar") {
    validate_separator_on_path(prob, feasible_paths, samples,
                               cptp::sep::MultistarSeparator{});
  }
  SECTION("Comb") {
    // Comb is enabled in path mode; this verifies validity if cuts are
    // produced.
    validate_separator_on_path(prob, feasible_paths, samples,
                               cptp::sep::CombSeparator{}, false);
  }
  SECTION("RGLM") {
    validate_separator_on_path(prob, feasible_paths, samples,
                               cptp::sep::RGLMSeparator{});
  }
}

TEST_CASE("Equivalent tour and s-t path trigger same cut families",
          "[separators][path][equiv]") {
  const auto inst = make_equivalent_instances();
  const auto samples = build_equivalent_samples(inst, 80);
  constexpr double tol = 1e-9;
  const double Q = inst.tour.capacity();
  double d_total = 0.0;
  for (int32_t i = 0; i < inst.tour.num_nodes(); ++i) {
    if (i != inst.tour.source()) d_total += inst.tour.demand(i);
  }

  // Customer subsets over {1,2,3,4}. In path model, 0=source and 5=target are
  // terminals.
  for (const auto& sample : samples) {
    for (int mask = 1; mask < (1 << 4); ++mask) {
      std::vector<uint8_t> in_tour(inst.tour.num_nodes(), 0);
      std::vector<uint8_t> in_path(inst.path.num_nodes(), 0);
      for (int b = 0; b < 4; ++b) {
        if (mask & (1 << b)) {
          const int32_t customer = b + 1;
          in_tour[customer] = 1;
          in_path[customer] = 1;
        }
      }

      const auto st_tour =
          compute_set_stats(inst.tour, sample.x_tour, sample.y_tour, in_tour);
      const auto st_path =
          compute_set_stats(inst.path, sample.x_path, sample.y_path, in_path);

      // Basic mapped quantities should match exactly.
      REQUIRE(std::abs(st_tour.dS - st_path.dS) <= tol);
      REQUIRE(std::abs(st_tour.served_demand - st_path.served_demand) <= tol);
      REQUIRE(std::abs(st_tour.cut_flow - st_path.cut_flow) <= tol);
      REQUIRE(std::abs(st_tour.star_sum - st_path.star_sum) <= tol);

      // SEC (for any target customer t in S): 2*y_t - x(delta(S))
      for (int b = 0; b < 4; ++b) {
        if (!(mask & (1 << b))) continue;
        const int32_t t = b + 1;
        const double sec_tour = 2.0 * sample.y_tour[t] - st_tour.cut_flow;
        const double sec_path = 2.0 * sample.y_path[t] - st_path.cut_flow;
        REQUIRE(std::abs(sec_tour - sec_path) <= tol);
      }

      // RCI violation
      const double Q_r = std::fmod(st_tour.dS, Q);
      if (Q_r > tol) {
        const double k = std::ceil(st_tour.dS / Q);
        if (k > 1.0) {
          const double rci_tour =
              2.0 * (k - st_tour.dS / Q_r) -
              (st_tour.cut_flow - (2.0 * st_tour.served_demand) / Q_r);
          const double rci_path =
              2.0 * (k - st_path.dS / Q_r) -
              (st_path.cut_flow - (2.0 * st_path.served_demand) / Q_r);
          REQUIRE(std::abs(rci_tour - rci_path) <= tol);
        }
      }

      // Multistar (GLM) LHS in <= form
      const double glm_tour = (2.0 / Q) * st_tour.served_demand -
                              (st_tour.cut_flow - (2.0 / Q) * st_tour.star_sum);
      const double glm_path = (2.0 / Q) * st_path.served_demand -
                              (st_path.cut_flow - (2.0 / Q) * st_path.star_sum);
      REQUIRE(std::abs(glm_tour - glm_path) <= tol);

      // RGLM violation
      const double alpha = 2.0 * d_total - st_tour.dS;
      const double r = std::fmod(alpha, Q);
      if (r > tol) {
        const double k = std::ceil(alpha / Q);
        if (k > 1.0) {
          const double beta_tour = (st_tour.dS - st_tour.served_demand) +
                                   2.0 * (d_total - st_tour.dS) -
                                   st_tour.star_sum;
          const double beta_path = (st_path.dS - st_path.served_demand) +
                                   2.0 * (d_total - st_path.dS) -
                                   st_path.star_sum;
          const double rglm_tour =
              (2.0 * k - 2.0 * beta_tour / r) - st_tour.cut_flow;
          const double rglm_path =
              (2.0 * k - 2.0 * beta_path / r) - st_path.cut_flow;
          REQUIRE(std::abs(rglm_tour - rglm_path) <= tol);
        }
      }
    }
  }
}

TEST_CASE("Equivalent tour and s-t path emit matching SPI customer cuts",
          "[spi][path][equiv]") {
  const auto inst = make_equivalent_instances();
  const auto samples = build_equivalent_samples(inst, 80);
  const auto ap_tour = compute_all_pairs(inst.tour);
  const auto ap_path = compute_all_pairs(inst.path);
  constexpr double ub = 90.0;
  constexpr double tol = 1e-6;

  auto normalize_spi = [](const std::vector<cptp::sep::Cut>& cuts,
                          int32_t y_offset) {
    std::set<uint32_t> keys;
    for (const auto& cut : cuts) {
      uint32_t mask = 0;
      for (size_t k = 0; k < cut.indices.size(); ++k) {
        REQUIRE(std::abs(cut.values[k] - 1.0) <= 1e-9);
        const int32_t node = cut.indices[k] - y_offset;
        REQUIRE(node >= 1);
        REQUIRE(node <= 4);
        mask |= (1u << static_cast<uint32_t>(node - 1));
      }
      const int rhs_int = static_cast<int>(std::lround(cut.rhs));
      REQUIRE(std::abs(cut.rhs - static_cast<double>(rhs_int)) <= 1e-6);
      keys.insert((static_cast<uint32_t>(rhs_int) << 8) | mask);
    }
    return keys;
  };

  for (const auto& sample : samples) {
    auto support_tour = build_support(inst.tour, sample.x_tour);
    auto support_path = build_support(inst.path, sample.x_path);

    cptp::sep::SeparationContext ctx_tour{
        .problem = inst.tour,
        .x_values = sample.x_tour,
        .y_values = sample.y_tour,
        .x_offset = 0,
        .y_offset = inst.tour.num_edges(),
        .tol = tol,
        .flow_tree = support_tour.tree.get(),
        .upper_bound = ub,
        .all_pairs = ap_tour,
    };
    cptp::sep::SeparationContext ctx_path{
        .problem = inst.path,
        .x_values = sample.x_path,
        .y_values = sample.y_path,
        .x_offset = 0,
        .y_offset = inst.path.num_edges(),
        .tol = tol,
        .flow_tree = support_path.tree.get(),
        .upper_bound = ub,
        .all_pairs = ap_path,
    };

    cptp::sep::SPISeparator spi_tour;
    cptp::sep::SPISeparator spi_path;
    const auto cuts_tour = spi_tour.separate(ctx_tour);
    const auto cuts_path = spi_path.separate(ctx_path);

    const auto keys_tour = normalize_spi(cuts_tour, inst.tour.num_edges());
    const auto keys_path = normalize_spi(cuts_path, inst.path.num_edges());
    REQUIRE(keys_tour == keys_path);
  }
}

TEST_CASE("Comb path emits valid cuts on crafted sample",
          "[comb][path][validity]") {
  const auto inst = make_comb_equivalent_instances();
  const auto sample = make_comb_equivalent_sample(inst);
  auto support = build_support(inst.path, sample.x_path);
  const int32_t m = inst.path.num_edges();
  const int32_t target = inst.path.target();
  const int32_t st_edge = edge_index(inst.path, inst.path.source(), target);
  REQUIRE(st_edge >= 0);
  cptp::sep::SeparationContext ctx{
      .problem = inst.path,
      .x_values = sample.x_path,
      .y_values = sample.y_path,
      .x_offset = 0,
      .y_offset = inst.path.num_edges(),
      .tol = 1e-6,
      .flow_tree = support.tree.get(),
  };
  cptp::sep::CombSeparator comb;
  const auto cuts = comb.separate(ctx);
  REQUIRE(!cuts.empty());

  const auto feasible_paths = enumerate_feasible_paths(inst.path);
  REQUIRE_FALSE(feasible_paths.empty());
  for (const auto& cut : cuts) {
    for (int32_t idx : cut.indices) {
      if (idx < m) {
        REQUIRE(idx != st_edge);
      } else {
        const int32_t node = idx - m;
        REQUIRE(node != target);
      }
    }
    for (const auto& p : feasible_paths) {
      const double lhs = evaluate_cut_lhs(cut, p.x_values, p.y_values);
      REQUIRE(lhs <= cut.rhs + 1e-6);
    }
  }
}

TEST_CASE(
    "Comb expression and separator behavior match across equivalent tour/path",
    "[comb][path][equiv]") {
  const auto inst = make_comb_equivalent_instances();
  const auto sample = make_comb_equivalent_sample(inst);

  // Expression-level equivalence on mapped handle/teeth.
  std::vector<uint8_t> H_tour(inst.tour.num_nodes(), 0);
  std::vector<uint8_t> H_path(inst.path.num_nodes(), 0);
  H_tour[0] = 1;
  H_path[0] = 1;
  H_path[inst.path.target()] = 1;
  for (int i : {1, 2, 3}) {
    H_tour[i] = 1;
    H_path[i] = 1;
  }
  std::vector<std::pair<int32_t, int32_t>> teeth = {{1, 4}, {2, 5}, {3, 6}};
  const double lhs_tour =
      compute_comb_lhs(inst.tour, sample.x_tour, sample.y_tour, H_tour, teeth);
  const double lhs_path = compute_comb_lhs_path_split(
      inst.path, sample.x_path, sample.y_path, H_path, teeth);
  const double rhs = 1.0;  // (t-1)/2 with t=3
  REQUIRE(std::abs(lhs_tour - lhs_path) <= 1e-9);
  REQUIRE(lhs_tour - rhs > 1e-6);

  // Separator-level check: both formulations should produce comb cuts, and
  // best violation should agree under the mapped LP sample.
  auto support_tour = build_support(inst.tour, sample.x_tour);
  auto support_path = build_support(inst.path, sample.x_path);
  cptp::sep::SeparationContext ctx_tour{
      .problem = inst.tour,
      .x_values = sample.x_tour,
      .y_values = sample.y_tour,
      .x_offset = 0,
      .y_offset = inst.tour.num_edges(),
      .tol = 1e-6,
      .flow_tree = support_tour.tree.get(),
  };
  cptp::sep::SeparationContext ctx_path{
      .problem = inst.path,
      .x_values = sample.x_path,
      .y_values = sample.y_path,
      .x_offset = 0,
      .y_offset = inst.path.num_edges(),
      .tol = 1e-6,
      .flow_tree = support_path.tree.get(),
  };
  cptp::sep::CombSeparator comb_tour;
  cptp::sep::CombSeparator comb_path;
  const auto cuts_tour = comb_tour.separate(ctx_tour);
  const auto cuts_path = comb_path.separate(ctx_path);
  REQUIRE(!cuts_tour.empty());
  REQUIRE(!cuts_path.empty());
  const int32_t m_path = inst.path.num_edges();
  const int32_t target = inst.path.target();
  const int32_t st_edge = edge_index(inst.path, inst.path.source(), target);
  REQUIRE(st_edge >= 0);
  for (const auto& cut : cuts_path) {
    for (int32_t idx : cut.indices) {
      if (idx < m_path) {
        REQUIRE(idx != st_edge);
      } else {
        const int32_t node = idx - m_path;
        REQUIRE(node != target);
      }
    }
  }

  auto max_violation = [](const std::vector<cptp::sep::Cut>& cuts) {
    double best = -std::numeric_limits<double>::infinity();
    for (const auto& c : cuts) best = std::max(best, c.violation);
    return best;
  };
  REQUIRE(std::abs(max_violation(cuts_tour) - max_violation(cuts_path)) <=
          1e-6);
}

TEST_CASE("Comb path ignores split-terminal artifact variables",
          "[comb][path][invariance]") {
  const auto inst = make_comb_equivalent_instances();
  const auto base = make_comb_equivalent_sample(inst);
  const int32_t m = inst.path.num_edges();
  const int32_t s = inst.path.source();
  const int32_t t = inst.path.target();
  const int32_t st_edge = edge_index(inst.path, s, t);
  REQUIRE(st_edge >= 0);

  auto run_comb_path = [&](const EquivalentSample& sample) {
    auto support = build_support(inst.path, sample.x_path);
    cptp::sep::SeparationContext ctx{
        .problem = inst.path,
        .x_values = sample.x_path,
        .y_values = sample.y_path,
        .x_offset = 0,
        .y_offset = m,
        .tol = 1e-6,
        .flow_tree = support.tree.get(),
    };
    cptp::sep::CombSeparator comb;
    return comb.separate(ctx);
  };

  const auto cuts_base = run_comb_path(base);
  REQUIRE(!cuts_base.empty());
  const auto sig_base = normalize_cut_signatures(cuts_base);

  // y_target is a split-terminal artifact and should not affect Comb path cuts.
  auto y_target_perturbed = base;
  y_target_perturbed.y_path[t] = 0.17;
  const auto sig_y_target =
      normalize_cut_signatures(run_comb_path(y_target_perturbed));
  REQUIRE(sig_y_target == sig_base);

  // x_{s,t} is a split-terminal artifact and should not affect Comb path cuts.
  auto st_perturbed = base;
  st_perturbed.x_path[st_edge] = 0.93;
  const auto sig_st = normalize_cut_signatures(run_comb_path(st_perturbed));
  REQUIRE(sig_st == sig_base);
}
