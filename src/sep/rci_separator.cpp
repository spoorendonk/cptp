#include "sep/rci_separator.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/gomory_hu.h"
#include "core/problem.h"

namespace cptp::sep {

namespace {

/// Working state for evaluating and refining an RCI candidate set.
struct RCICandidate {
  std::vector<bool> in_S;  // membership flags
  double d_S;              // total demand of S
  double served_demand;    // sum q_i * y_i for i in S
  double cut_flow;         // x(delta(S))
  double inside_flow;      // x(E(S))
};

/// Compute violation of the CPTP RCI for the given candidate.
/// Returns violation (positive = violated).
double compute_violation(const RCICandidate& c, double Q, double tol) {
  double Q_r = std::fmod(c.d_S, Q);
  if (Q_r <= tol) return -1.0;

  double k = std::ceil(c.d_S / Q);
  if (k <= 1.0) return -1.0;

  double rhs = 2.0 * (k - c.d_S / Q_r);
  double lhs = c.cut_flow - (2.0 * c.served_demand) / Q_r;
  return rhs - lhs;
}

/// Initialize a candidate from a GH tree cut partition.
/// For non-bottleneck cuts, the depot may be on either side.
/// S is always the side NOT containing the depot.
RCICandidate make_candidate(const gomory_hu_tree& tree, int32_t cut_node,
                            const Problem& prob, const SeparationContext& ctx) {
  const auto& graph = prob.graph();
  const int32_t n = prob.num_nodes();
  const int32_t depot = prob.source();
  const bool is_tour = prob.is_tour();
  const int32_t path_target = prob.target();

  bool depot_on_source = !tree.on_root_side(cut_node, depot);

  RCICandidate c;
  c.in_S.resize(n, false);
  c.d_S = 0.0;
  c.served_demand = 0.0;
  c.cut_flow = 0.0;
  c.inside_flow = 0.0;

  for (int32_t i = 0; i < n; ++i) {
    bool on_source_side = !tree.on_root_side(cut_node, i);
    bool in_S = (depot_on_source) ? !on_source_side : on_source_side;

    if (in_S && i != depot && (is_tour || i != path_target)) {
      c.in_S[i] = true;
      c.d_S += prob.demand(i);
      c.served_demand += prob.demand(i) * ctx.y_values[i];
    }
  }

  for (auto e : graph.edges()) {
    int32_t u = graph.edge_source(e);
    int32_t v = graph.edge_target(e);
    bool u_in = c.in_S[u];
    bool v_in = c.in_S[v];
    double x_e = ctx.x_values[e];
    if (u_in && v_in) {
      c.inside_flow += x_e;
    } else if (u_in != v_in) {
      c.cut_flow += x_e;
    }
  }

  return c;
}

/// Incrementally update candidate when adding node j to S.
void add_node(RCICandidate& c, int32_t j, const Problem& prob,
              const SeparationContext& ctx) {
  const auto& graph = prob.graph();
  c.in_S[j] = true;
  c.d_S += prob.demand(j);
  c.served_demand += prob.demand(j) * ctx.y_values[j];

  for (int32_t e : graph.incident_edges(j)) {
    int32_t other = graph.other_endpoint(e, j);
    double x_e = ctx.x_values[e];
    if (c.in_S[other]) {
      c.cut_flow -= x_e;
      c.inside_flow += x_e;
    } else {
      c.cut_flow += x_e;
    }
  }
}

/// Incrementally update candidate when removing node j from S.
void remove_node(RCICandidate& c, int32_t j, const Problem& prob,
                 const SeparationContext& ctx) {
  const auto& graph = prob.graph();
  c.in_S[j] = false;
  c.d_S -= prob.demand(j);
  c.served_demand -= prob.demand(j) * ctx.y_values[j];

  for (int32_t e : graph.incident_edges(j)) {
    int32_t other = graph.other_endpoint(e, j);
    double x_e = ctx.x_values[e];
    if (c.in_S[other]) {
      c.inside_flow -= x_e;
      c.cut_flow += x_e;
    } else {
      c.cut_flow -= x_e;
    }
  }
}

/// Run add/drop local search on a candidate.
/// Returns the best violation found.
double add_drop_search(RCICandidate& c, const Problem& prob,
                       const SeparationContext& ctx, double Q, double tol) {
  const auto& graph = prob.graph();
  const int32_t n = prob.num_nodes();
  const int32_t depot = prob.source();
  const bool is_tour = prob.is_tour();
  const int32_t path_target = prob.target();

  double best_viol = compute_violation(c, Q, tol);
  bool improved = true;

  while (improved) {
    improved = false;

    // --- Drop phase ---
    bool drop_improved = true;
    while (drop_improved) {
      drop_improved = false;
      int32_t best_drop = -1;
      double best_drop_viol = best_viol;

      for (int32_t i = 0; i < n; ++i) {
        if (!c.in_S[i] || i == depot || (!is_tour && i == path_target))
          continue;

        remove_node(c, i, prob, ctx);
        double v = compute_violation(c, Q, tol);
        if (v > best_drop_viol + tol) {
          best_drop_viol = v;
          best_drop = i;
        }
        add_node(c, i, prob, ctx);
      }

      if (best_drop >= 0) {
        remove_node(c, best_drop, prob, ctx);
        best_viol = best_drop_viol;
        drop_improved = true;
        improved = true;
      }
    }

    // --- Add phase ---
    bool add_improved = true;
    while (add_improved) {
      add_improved = false;
      int32_t best_add = -1;
      double best_add_viol = best_viol;

      for (int32_t i = 0; i < n; ++i) {
        if (c.in_S[i] || i == depot || (!is_tour && i == path_target)) continue;

        bool adjacent = false;
        for (int32_t e : graph.incident_edges(i)) {
          int32_t other = graph.other_endpoint(e, i);
          if (c.in_S[other] && ctx.x_values[e] > tol) {
            adjacent = true;
            break;
          }
        }
        if (!adjacent) continue;

        add_node(c, i, prob, ctx);
        double v = compute_violation(c, Q, tol);
        if (v > best_add_viol + tol) {
          best_add_viol = v;
          best_add = i;
        }
        remove_node(c, i, prob, ctx);
      }

      if (best_add >= 0) {
        add_node(c, best_add, prob, ctx);
        best_viol = best_add_viol;
        add_improved = true;
        improved = true;
      }
    }
  }

  return best_viol;
}

/// Build a Cut from a finalized candidate, choosing the sparser form.
Cut build_cut(const RCICandidate& c, const Problem& prob,
              const SeparationContext& ctx, double Q, double tol) {
  const auto& graph = prob.graph();
  const int32_t n = prob.num_nodes();

  double Q_r = std::fmod(c.d_S, Q);
  double k = std::ceil(c.d_S / Q);

  // Count |E(S)| and |δ(S)| to choose sparser form
  int32_t n_inside = 0;
  int32_t n_delta = 0;
  for (auto e : graph.edges()) {
    int32_t u = graph.edge_source(e);
    int32_t v = graph.edge_target(e);
    bool u_in = c.in_S[u];
    bool v_in = c.in_S[v];
    if (u_in && v_in)
      n_inside++;
    else if (u_in != v_in)
      n_delta++;
  }

  bool use_inside_form = (n_inside < n_delta);

  Cut cut;
  cut.violation = compute_violation(c, Q, tol);

  if (use_inside_form) {
    // Inside form (via degree substitution x(δ({i}))=2y_i):
    //   x(E(S)) - Σ_{i∈S} (1 - q_i/Q_r)·y_i  ≤  d(S)/Q_r - k
    for (auto e : graph.edges()) {
      int32_t u = graph.edge_source(e);
      int32_t v = graph.edge_target(e);
      if (c.in_S[u] && c.in_S[v]) {
        cut.indices.push_back(ctx.x_offset + e);
        cut.values.push_back(1.0);
      }
    }

    for (int32_t i = 0; i < n; ++i) {
      if (!c.in_S[i]) continue;
      double q_i = prob.demand(i);
      double coeff = -(1.0 - q_i / Q_r);
      if (std::abs(coeff) > tol) {
        cut.indices.push_back(ctx.y_offset + i);
        cut.values.push_back(coeff);
      }
    }

    cut.rhs = c.d_S / Q_r - k;
  } else {
    // Cut form: (2/Q_r)·Σ q_i·y_i - x(δ(S)) ≤ -2·(k - d(S)/Q_r)
    for (auto e : graph.edges()) {
      int32_t u = graph.edge_source(e);
      int32_t v = graph.edge_target(e);
      if (c.in_S[u] != c.in_S[v]) {
        cut.indices.push_back(ctx.x_offset + e);
        cut.values.push_back(-1.0);
      }
    }

    for (int32_t i = 0; i < n; ++i) {
      if (!c.in_S[i]) continue;
      double q_i = prob.demand(i);
      if (q_i <= 0) continue;

      double coeff = 2.0 * q_i / Q_r;
      cut.indices.push_back(ctx.y_offset + i);
      cut.values.push_back(coeff);
    }

    cut.rhs = -2.0 * (k - c.d_S / Q_r);
  }

  return cut;
}

}  // namespace

std::vector<Cut> RCISeparator::separate(const SeparationContext& ctx) {
  const auto& prob = ctx.problem;
  const int32_t n = prob.num_nodes();
  const double tol = ctx.tol;
  const double Q = prob.capacity();
  const int32_t depot = prob.source();
  const bool is_tour = prob.is_tour();
  const int32_t path_target = prob.target();

  if (Q <= 0 || !ctx.flow_tree) return {};

  const auto& tree = *ctx.flow_tree;

  // Step 1: Collect all unique GH tree cuts across all target paths.
  std::vector<bool> evaluated(n, false);

  std::vector<gomory_hu_tree::CutRef> candidates;
  for (int32_t target = 0; target < n; ++target) {
    if (target == depot || (!is_tour && target == path_target) ||
        ctx.y_values[target] <= tol) {
      continue;
    }

    auto path_cuts = tree.all_cuts_on_path(target);
    for (auto& cr : path_cuts) {
      if (!evaluated[cr.node]) {
        evaluated[cr.node] = true;
        candidates.push_back(cr);
      }
    }
  }

  if (candidates.empty()) return {};

  // Step 2: Parallel evaluation and add/drop local search.
  // Each candidate produces at most one cut.
  struct CutResult {
    Cut cut;
    double violation;
  };
  std::vector<std::vector<CutResult>> per_candidate(candidates.size());

  tbb::parallel_for(tbb::blocked_range<size_t>(0, candidates.size()),
                    [&](const tbb::blocked_range<size_t>& range) {
                      for (size_t idx = range.begin(); idx < range.end();
                           ++idx) {
                        auto& cr = candidates[idx];

                        auto c = make_candidate(tree, cr.node, prob, ctx);

                        double Q_r = std::fmod(c.d_S, Q);
                        if (Q_r <= tol) continue;
                        double k = std::ceil(c.d_S / Q);
                        if (k <= 1.0) continue;

                        double viol = compute_violation(c, Q, tol);

                        // Run add/drop on violated candidates to strengthen.
                        if (viol > tol) {
                          viol = add_drop_search(c, prob, ctx, Q, tol);
                        }

                        if (viol > tol) {
                          auto cut = build_cut(c, prob, ctx, Q, tol);

                          per_candidate[idx].push_back({std::move(cut), viol});
                        }
                      }
                    });

  // Step 3: Collect, sort by violation, and cap output.
  std::vector<CutResult> results;
  for (auto& pc : per_candidate) {
    for (auto& cr : pc) {
      results.push_back(std::move(cr));
    }
  }

  // Sort by violation descending — most violated first.
  // stable_sort preserves deterministic input order for equal violations.
  std::stable_sort(results.begin(), results.end(),
                   [](const CutResult& a, const CutResult& b) {
                     return a.violation > b.violation;
                   });

  // Cap output: keep at most max_cuts per round. This avoids flooding
  // the LP with weak cuts that slow down the simplex.
  const size_t max_cuts = std::max(5, n / 5);

  std::vector<Cut> all_cuts;
  all_cuts.reserve(std::min(results.size(), max_cuts));
  for (size_t i = 0; i < std::min(results.size(), max_cuts); ++i) {
    all_cuts.push_back(std::move(results[i].cut));
  }

  return all_cuts;
}

}  // namespace cptp::sep
