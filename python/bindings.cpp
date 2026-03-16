#include <memory>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "core/io.h"
#include "core/problem.h"
#include "core/solution.h"
#include "model/model.h"

namespace nb = nanobind;
using namespace nb::literals;

// Helper: wrap a const std::vector<T>& as a read-only numpy view (no copy).
// The parent python object must stay alive (use rv_policy::reference_internal).
template <typename T>
static nb::ndarray<nb::numpy, const T, nb::shape<-1>> vec_view_numpy(
    const std::vector<T>& v, nb::handle parent) {
  return nb::ndarray<nb::numpy, const T, nb::shape<-1>>(
      const_cast<T*>(v.data()), {v.size()}, parent);
}

NB_MODULE(_cptp, m) {
  m.doc() = "CPTP Branch-and-Cut Solver";

  // --- Status enum ---
  nb::enum_<cptp::SolveResult::Status>(m, "Status")
      .value("Optimal", cptp::SolveResult::Status::Optimal)
      .value("Feasible", cptp::SolveResult::Status::Feasible)
      .value("Infeasible", cptp::SolveResult::Status::Infeasible)
      .value("Unbounded", cptp::SolveResult::Status::Unbounded)
      .value("TimeLimit", cptp::SolveResult::Status::TimeLimit)
      .value("Error", cptp::SolveResult::Status::Error);

  // --- SeparatorStats ---
  nb::class_<cptp::SeparatorStats>(m, "SeparatorStats")
      .def_ro("cuts_added", &cptp::SeparatorStats::cuts_added)
      .def_ro("rounds_called", &cptp::SeparatorStats::rounds_called)
      .def_ro("time_seconds", &cptp::SeparatorStats::time_seconds);

  // --- SolveResult ---
  // tour and tour_arcs are returned as numpy arrays (zero-copy move from C++
  // vector).
  nb::class_<cptp::SolveResult>(m, "SolveResult")
      .def_ro("status", &cptp::SolveResult::status)
      .def_ro("objective", &cptp::SolveResult::objective)
      .def_ro("bound", &cptp::SolveResult::bound)
      .def_ro("gap", &cptp::SolveResult::gap)
      .def_ro("time_seconds", &cptp::SolveResult::time_seconds)
      .def_ro("nodes", &cptp::SolveResult::nodes)
      .def_prop_ro("tour",
                   [](cptp::SolveResult& self) {
                     return vec_view_numpy<int32_t>(self.tour, nb::find(self));
                   })
      .def_prop_ro("tour_arcs",
                   [](cptp::SolveResult& self) {
                     return vec_view_numpy<int32_t>(self.tour_arcs,
                                                    nb::find(self));
                   })
      .def_ro("total_cuts", &cptp::SolveResult::total_cuts)
      .def_ro("separation_rounds", &cptp::SolveResult::separation_rounds)
      .def_ro("separator_stats", &cptp::SolveResult::separator_stats)
      .def("is_optimal", &cptp::SolveResult::is_optimal)
      .def("has_solution", &cptp::SolveResult::has_solution);

  // --- Problem ---
  // Vector accessors return numpy views (zero-copy, backed by the Problem
  // object).
  nb::class_<cptp::Problem>(m, "Problem")
      .def(nb::init<>())
      .def_rw("name", &cptp::Problem::name)
      .def_prop_ro("num_nodes", &cptp::Problem::num_nodes)
      .def_prop_ro("num_edges", &cptp::Problem::num_edges)
      .def_prop_ro("source", &cptp::Problem::source)
      .def_prop_ro("target", &cptp::Problem::target)
      .def_prop_ro("capacity", &cptp::Problem::capacity)
      .def_prop_ro("is_tour", &cptp::Problem::is_tour)
      .def_prop_ro("edge_costs",
                   [](cptp::Problem& self) {
                     return vec_view_numpy<double>(self.edge_costs(),
                                                   nb::find(self));
                   })
      .def_prop_ro("profits",
                   [](cptp::Problem& self) {
                     return vec_view_numpy<double>(self.profits(),
                                                   nb::find(self));
                   })
      .def_prop_ro("demands",
                   [](cptp::Problem& self) {
                     return vec_view_numpy<double>(self.demands(),
                                                   nb::find(self));
                   })
      .def(
          "graph_edges",
          [](const cptp::Problem& self) {
            const auto& g = self.graph();
            size_t m = static_cast<size_t>(self.num_edges());
            auto data = std::make_unique<std::vector<int32_t>>(m * 2);
            size_t idx = 0;
            for (auto e : g.edges()) {
              (*data)[idx++] = g.edge_source(e);
              (*data)[idx++] = g.edge_target(e);
            }
            auto* raw = data.release();
            nb::capsule owner(raw, [](void* p) noexcept {
              delete static_cast<std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t, nb::shape<-1, 2>>(
                raw->data(), {m, 2}, std::move(owner));
          },
          "Return (m, 2) numpy array of (tail, head) pairs")
      .def(
          "__init__",
          [](cptp::Problem* self, int32_t num_nodes,
             nb::ndarray<int32_t, nb::shape<-1, 2>> edges,
             nb::ndarray<double, nb::shape<-1>> edge_costs,
             nb::ndarray<double, nb::shape<-1>> profits,
             nb::ndarray<double, nb::shape<-1>> demands, double capacity,
             int32_t source, int32_t target, const std::string& name) {
            new (self) cptp::Problem();
            auto e_view = edges.view();
            int32_t m = static_cast<int32_t>(edges.shape(0));
            std::vector<cptp::Edge> edge_vec(m);
            for (int32_t i = 0; i < m; ++i)
              edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            auto ec_view = edge_costs.view();
            auto p_view = profits.view();
            auto d_view = demands.view();
            self->build(
                num_nodes, edge_vec,
                {ec_view.data(), static_cast<size_t>(edge_costs.shape(0))},
                {p_view.data(), static_cast<size_t>(profits.shape(0))},
                {d_view.data(), static_cast<size_t>(demands.shape(0))},
                capacity, source, target);
            self->name = name;
          },
          "num_nodes"_a, "edges"_a, "edge_costs"_a, "profits"_a, "demands"_a,
          "capacity"_a, "source"_a = 0, "target"_a = 0, "name"_a = "");

  // --- Model ---
  // Input arrays are read directly from numpy (zero-copy read via span).
  nb::class_<cptp::Model>(m, "Model")
      .def(nb::init<>())
      .def("set_problem", &cptp::Model::set_problem, "problem"_a)
      .def("problem", &cptp::Model::problem, nb::rv_policy::reference_internal)
      .def(
          "set_graph",
          [](cptp::Model& self, int32_t n,
             nb::ndarray<int32_t, nb::shape<-1, 2>> edges,
             nb::ndarray<double, nb::shape<-1>> costs) {
            auto e_view = edges.view();
            auto c_view = costs.view();
            int32_t num_edges = static_cast<int32_t>(edges.shape(0));

            std::vector<cptp::Edge> edge_vec(num_edges);
            for (int32_t i = 0; i < num_edges; ++i) {
              edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            }

            self.set_graph(n, edge_vec,
                           {c_view.data(), static_cast<size_t>(num_edges)});
          },
          "num_nodes"_a, "edges"_a, "edge_costs"_a)
      .def("set_source", &cptp::Model::set_source)
      .def("set_target", &cptp::Model::set_target)
      .def("set_depot", &cptp::Model::set_depot)
      .def(
          "set_profits",
          [](cptp::Model& self, nb::ndarray<double, nb::shape<-1>> profits) {
            auto view = profits.view();
            self.set_profits(
                {view.data(), static_cast<size_t>(profits.shape(0))});
          },
          "profits"_a)
      .def(
          "add_capacity_resource",
          [](cptp::Model& self, nb::ndarray<double, nb::shape<-1>> demands,
             double limit) {
            auto view = demands.view();
            self.add_capacity_resource(
                {view.data(), static_cast<size_t>(demands.shape(0))}, limit);
          },
          "demands"_a, "limit"_a)
      .def(
          "solve",
          [](cptp::Model& self, const cptp::SolverOptions& options) {
            return self.solve(options);
          },
          "options"_a = cptp::SolverOptions{});

  m.attr("has_highs") = true;

  // IO functions (always available)
  m.def(
      "load", [](const std::string& path) { return cptp::io::load(path); },
      "path"_a, "Load an CPTP instance from file (auto-detect format)");
}
