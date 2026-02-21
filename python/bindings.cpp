#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "model/model.h"
#include "core/io.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_rcspp_bac, m) {
    m.doc() = "RCSPP Branch-and-Cut Solver";

    // --- Status enum ---
    nb::enum_<rcspp::SolveResult::Status>(m, "Status")
        .value("Optimal", rcspp::SolveResult::Status::Optimal)
        .value("Feasible", rcspp::SolveResult::Status::Feasible)
        .value("Infeasible", rcspp::SolveResult::Status::Infeasible)
        .value("Unbounded", rcspp::SolveResult::Status::Unbounded)
        .value("TimeLimit", rcspp::SolveResult::Status::TimeLimit)
        .value("Error", rcspp::SolveResult::Status::Error);

    // --- SeparatorStats ---
    nb::class_<rcspp::SeparatorStats>(m, "SeparatorStats")
        .def_ro("cuts_added", &rcspp::SeparatorStats::cuts_added)
        .def_ro("rounds_called", &rcspp::SeparatorStats::rounds_called)
        .def_ro("time_seconds", &rcspp::SeparatorStats::time_seconds);

    // --- SolveResult ---
    nb::class_<rcspp::SolveResult>(m, "SolveResult")
        .def_ro("status", &rcspp::SolveResult::status)
        .def_ro("objective", &rcspp::SolveResult::objective)
        .def_ro("bound", &rcspp::SolveResult::bound)
        .def_ro("gap", &rcspp::SolveResult::gap)
        .def_ro("time_seconds", &rcspp::SolveResult::time_seconds)
        .def_ro("nodes", &rcspp::SolveResult::nodes)
        .def_ro("tour", &rcspp::SolveResult::tour)
        .def_ro("tour_arcs", &rcspp::SolveResult::tour_arcs)
        .def_ro("total_cuts", &rcspp::SolveResult::total_cuts)
        .def_ro("separation_rounds", &rcspp::SolveResult::separation_rounds)
        .def_ro("separator_stats", &rcspp::SolveResult::separator_stats)
        .def("is_optimal", &rcspp::SolveResult::is_optimal)
        .def("has_solution", &rcspp::SolveResult::has_solution);

    // --- Problem ---
    nb::class_<rcspp::Problem>(m, "Problem")
        .def(nb::init<>())
        .def_rw("name", &rcspp::Problem::name)
        .def_prop_ro("num_nodes", &rcspp::Problem::num_nodes)
        .def_prop_ro("num_edges", &rcspp::Problem::num_edges)
        .def_prop_ro("source", &rcspp::Problem::source)
        .def_prop_ro("target", &rcspp::Problem::target)
        .def_prop_ro("capacity", &rcspp::Problem::capacity)
        .def_prop_ro("is_tour", &rcspp::Problem::is_tour)
        .def_prop_ro("edge_costs", &rcspp::Problem::edge_costs)
        .def_prop_ro("profits", &rcspp::Problem::profits)
        .def_prop_ro("demands", &rcspp::Problem::demands)
        .def("graph_edges", [](const rcspp::Problem& self) {
            std::vector<std::pair<int32_t, int32_t>> edges;
            const auto& g = self.graph();
            for (auto e : g.edges())
                edges.push_back({g.edge_source(e), g.edge_target(e)});
            return edges;
        }, "Return list of (tail, head) pairs for all edges")
        .def("__init__", [](rcspp::Problem* self,
                            int32_t num_nodes,
                            const std::vector<std::pair<int32_t, int32_t>>& edges,
                            const std::vector<double>& edge_costs,
                            const std::vector<double>& profits,
                            const std::vector<double>& demands,
                            double capacity,
                            int32_t source,
                            int32_t target,
                            const std::string& name) {
            new (self) rcspp::Problem();
            std::vector<rcspp::Edge> edge_vec(edges.size());
            for (size_t i = 0; i < edges.size(); ++i)
                edge_vec[i] = {edges[i].first, edges[i].second};
            self->build(num_nodes, edge_vec, edge_costs, profits, demands, capacity, source, target);
            self->name = name;
        }, "num_nodes"_a, "edges"_a, "edge_costs"_a, "profits"_a,
           "demands"_a, "capacity"_a, "source"_a = 0, "target"_a = 0,
           "name"_a = "");

    // --- Model ---
    nb::class_<rcspp::Model>(m, "Model")
        .def(nb::init<>())
        .def("set_problem", &rcspp::Model::set_problem, "problem"_a)
        .def("problem", &rcspp::Model::problem, nb::rv_policy::reference_internal)
        .def("set_graph", [](rcspp::Model& self, int32_t n,
                              nb::ndarray<int32_t, nb::shape<nb::any, 2>> edges,
                              nb::ndarray<double, nb::shape<nb::any>> costs) {
            auto e_view = edges.view();
            auto c_view = costs.view();
            int32_t num_edges = static_cast<int32_t>(edges.shape(0));

            std::vector<rcspp::Edge> edge_vec(num_edges);
            for (int32_t i = 0; i < num_edges; ++i) {
                edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            }

            std::vector<double> cost_vec(c_view.data(), c_view.data() + num_edges);
            self.set_graph(n, edge_vec, cost_vec);
        }, "num_nodes"_a, "edges"_a, "edge_costs"_a)
        .def("set_source", &rcspp::Model::set_source)
        .def("set_target", &rcspp::Model::set_target)
        .def("set_depot", &rcspp::Model::set_depot)
        .def("set_profits", [](rcspp::Model& self,
                                nb::ndarray<double, nb::shape<nb::any>> profits) {
            auto view = profits.view();
            std::vector<double> vec(view.data(), view.data() + profits.shape(0));
            self.set_profits(vec);
        }, "profits"_a)
        .def("add_capacity_resource", [](rcspp::Model& self,
                                          nb::ndarray<double, nb::shape<nb::any>> demands,
                                          double limit) {
            auto view = demands.view();
            std::vector<double> vec(view.data(), view.data() + demands.shape(0));
            self.add_capacity_resource(vec, limit);
        }, "demands"_a, "limit"_a)
        .def("solve", [](rcspp::Model& self,
                          const rcspp::SolverOptions& options) {
            return self.solve(options);
        }, "options"_a = rcspp::SolverOptions{});

    // IO functions
    m.def("load", [](const std::string& path) {
        return rcspp::io::load(path);
    }, "path"_a, "Load an RCSPP instance from file (auto-detect format)");
}
