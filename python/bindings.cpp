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

NB_MODULE(_cptp, m) {
    m.doc() = "CPTP Branch-and-Cut Solver";

    nb::class_<cptp::SolveResult>(m, "SolveResult")
        .def_ro("objective", &cptp::SolveResult::objective)
        .def_ro("bound", &cptp::SolveResult::bound)
        .def_ro("gap", &cptp::SolveResult::gap)
        .def_ro("time_seconds", &cptp::SolveResult::time_seconds)
        .def_ro("nodes", &cptp::SolveResult::nodes)
        .def_ro("tour", &cptp::SolveResult::tour)
        .def_ro("tour_arcs", &cptp::SolveResult::tour_arcs)
        .def("is_optimal", &cptp::SolveResult::is_optimal)
        .def("has_solution", &cptp::SolveResult::has_solution);

    nb::class_<cptp::Model>(m, "Model")
        .def(nb::init<>())
        .def("set_graph", [](cptp::Model& self, int32_t n,
                              nb::ndarray<int32_t, nb::shape<nb::any, 2>> edges,
                              nb::ndarray<double, nb::shape<nb::any>> costs) {
            auto e_view = edges.view();
            auto c_view = costs.view();
            int32_t num_edges = static_cast<int32_t>(edges.shape(0));

            std::vector<cptp::Edge> edge_vec(num_edges);
            for (int32_t i = 0; i < num_edges; ++i) {
                edge_vec[i] = {e_view(i, 0), e_view(i, 1)};
            }

            std::vector<double> cost_vec(c_view.data(), c_view.data() + num_edges);
            self.set_graph(n, edge_vec, cost_vec);
        }, "num_nodes"_a, "edges"_a, "edge_costs"_a)
        .def("set_source", &cptp::Model::set_source)
        .def("set_target", &cptp::Model::set_target)
        .def("set_depot", &cptp::Model::set_depot)
        .def("set_profits", [](cptp::Model& self,
                                nb::ndarray<double, nb::shape<nb::any>> profits) {
            auto view = profits.view();
            std::vector<double> vec(view.data(), view.data() + profits.shape(0));
            self.set_profits(vec);
        }, "profits"_a)
        .def("add_capacity_resource", [](cptp::Model& self,
                                          nb::ndarray<double, nb::shape<nb::any>> demands,
                                          double limit) {
            auto view = demands.view();
            std::vector<double> vec(view.data(), view.data() + demands.shape(0));
            self.add_capacity_resource(vec, limit);
        }, "demands"_a, "limit"_a)
        .def("solve", [](cptp::Model& self,
                          const cptp::SolverOptions& options) {
            return self.solve(options);
        }, "options"_a = cptp::SolverOptions{});

    // IO functions
    m.def("load", [](const std::string& path) {
        return cptp::io::load(path);
    }, "path"_a, "Load a CPTP instance from file (auto-detect format)");
}
