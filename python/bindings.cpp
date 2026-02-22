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

// Helper: wrap a const std::vector<T>& as a read-only numpy view (no copy).
// The parent python object must stay alive (use rv_policy::reference_internal).
template <typename T>
static nb::ndarray<nb::numpy, const T, nb::shape<-1>>
vec_view_numpy(const std::vector<T>& v, nb::handle parent) {
    return nb::ndarray<nb::numpy, const T, nb::shape<-1>>(
        const_cast<T*>(v.data()), {v.size()}, parent);
}

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
    // tour and tour_arcs are returned as numpy arrays (zero-copy move from C++ vector).
    nb::class_<rcspp::SolveResult>(m, "SolveResult")
        .def_ro("status", &rcspp::SolveResult::status)
        .def_ro("objective", &rcspp::SolveResult::objective)
        .def_ro("bound", &rcspp::SolveResult::bound)
        .def_ro("gap", &rcspp::SolveResult::gap)
        .def_ro("time_seconds", &rcspp::SolveResult::time_seconds)
        .def_ro("nodes", &rcspp::SolveResult::nodes)
        .def_prop_ro("tour", [](rcspp::SolveResult& self) {
            return vec_view_numpy<int32_t>(self.tour, nb::find(self));
        })
        .def_prop_ro("tour_arcs", [](rcspp::SolveResult& self) {
            return vec_view_numpy<int32_t>(self.tour_arcs, nb::find(self));
        })
        .def_ro("total_cuts", &rcspp::SolveResult::total_cuts)
        .def_ro("separation_rounds", &rcspp::SolveResult::separation_rounds)
        .def_ro("separator_stats", &rcspp::SolveResult::separator_stats)
        .def("is_optimal", &rcspp::SolveResult::is_optimal)
        .def("has_solution", &rcspp::SolveResult::has_solution);

    // --- Problem ---
    // Vector accessors return numpy views (zero-copy, backed by the Problem object).
    nb::class_<rcspp::Problem>(m, "Problem")
        .def(nb::init<>())
        .def_rw("name", &rcspp::Problem::name)
        .def_prop_ro("num_nodes", &rcspp::Problem::num_nodes)
        .def_prop_ro("num_edges", &rcspp::Problem::num_edges)
        .def_prop_ro("source", &rcspp::Problem::source)
        .def_prop_ro("target", &rcspp::Problem::target)
        .def_prop_ro("capacity", &rcspp::Problem::capacity)
        .def_prop_ro("is_tour", &rcspp::Problem::is_tour)
        .def_prop_ro("edge_costs", [](rcspp::Problem& self) {
            return vec_view_numpy<double>(self.edge_costs(), nb::find(self));
        })
        .def_prop_ro("profits", [](rcspp::Problem& self) {
            return vec_view_numpy<double>(self.profits(), nb::find(self));
        })
        .def_prop_ro("demands", [](rcspp::Problem& self) {
            return vec_view_numpy<double>(self.demands(), nb::find(self));
        })
        .def("graph_edges", [](const rcspp::Problem& self) {
            const auto& g = self.graph();
            size_t m = static_cast<size_t>(self.num_edges());
            auto* data = new std::vector<int32_t>(m * 2);
            size_t idx = 0;
            for (auto e : g.edges()) {
                (*data)[idx++] = g.edge_source(e);
                (*data)[idx++] = g.edge_target(e);
            }
            nb::capsule owner(data, [](void* p) noexcept {
                delete static_cast<std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t, nb::shape<-1, 2>>(
                data->data(), {m, 2}, std::move(owner));
        }, "Return (m, 2) numpy array of (tail, head) pairs")
        .def("__init__", [](rcspp::Problem* self,
                            int32_t num_nodes,
                            nb::ndarray<int32_t, nb::shape<-1, 2>> edges,
                            nb::ndarray<double, nb::shape<-1>> edge_costs,
                            nb::ndarray<double, nb::shape<-1>> profits,
                            nb::ndarray<double, nb::shape<-1>> demands,
                            double capacity,
                            int32_t source,
                            int32_t target,
                            const std::string& name) {
            new (self) rcspp::Problem();
            auto e_view = edges.view();
            int32_t m = static_cast<int32_t>(edges.shape(0));
            std::vector<rcspp::Edge> edge_vec(m);
            for (int32_t i = 0; i < m; ++i)
                edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            auto ec_view = edge_costs.view();
            auto p_view = profits.view();
            auto d_view = demands.view();
            self->build(num_nodes, edge_vec,
                        {ec_view.data(), static_cast<size_t>(edge_costs.shape(0))},
                        {p_view.data(), static_cast<size_t>(profits.shape(0))},
                        {d_view.data(), static_cast<size_t>(demands.shape(0))},
                        capacity, source, target);
            self->name = name;
        }, "num_nodes"_a, "edges"_a, "edge_costs"_a, "profits"_a,
           "demands"_a, "capacity"_a, "source"_a = 0, "target"_a = 0,
           "name"_a = "");

    // --- Model ---
    // Input arrays are read directly from numpy (zero-copy read via span).
    nb::class_<rcspp::Model>(m, "Model")
        .def(nb::init<>())
        .def("set_problem", &rcspp::Model::set_problem, "problem"_a)
        .def("problem", &rcspp::Model::problem, nb::rv_policy::reference_internal)
        .def("set_graph", [](rcspp::Model& self, int32_t n,
                              nb::ndarray<int32_t, nb::shape<-1, 2>> edges,
                              nb::ndarray<double, nb::shape<-1>> costs) {
            auto e_view = edges.view();
            auto c_view = costs.view();
            int32_t num_edges = static_cast<int32_t>(edges.shape(0));

            std::vector<rcspp::Edge> edge_vec(num_edges);
            for (int32_t i = 0; i < num_edges; ++i) {
                edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            }

            self.set_graph(n, edge_vec,
                           {c_view.data(), static_cast<size_t>(num_edges)});
        }, "num_nodes"_a, "edges"_a, "edge_costs"_a)
        .def("set_source", &rcspp::Model::set_source)
        .def("set_target", &rcspp::Model::set_target)
        .def("set_depot", &rcspp::Model::set_depot)
        .def("set_profits", [](rcspp::Model& self,
                                nb::ndarray<double, nb::shape<-1>> profits) {
            auto view = profits.view();
            self.set_profits({view.data(), static_cast<size_t>(profits.shape(0))});
        }, "profits"_a)
        .def("add_capacity_resource", [](rcspp::Model& self,
                                          nb::ndarray<double, nb::shape<-1>> demands,
                                          double limit) {
            auto view = demands.view();
            self.add_capacity_resource(
                {view.data(), static_cast<size_t>(demands.shape(0))}, limit);
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
