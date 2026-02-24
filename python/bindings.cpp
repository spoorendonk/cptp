#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "core/io.h"
#include "core/problem.h"
#include "core/solution.h"
#include "heuristic/primal_heuristic.h"
#include "preprocess/bound_propagator.h"
#include "preprocess/edge_elimination.h"
#include "sep/comb_separator.h"
#include "sep/cut.h"
#include "sep/multistar_separator.h"
#include "sep/rci_separator.h"
#include "sep/rglm_separator.h"
#include "sep/sec_separator.h"
#include "sep/separation_oracle.h"
#include "sep/separator.h"

#ifdef RCSPP_HAS_HIGHS
#include "model/model.h"
#endif

namespace nb = nanobind;
using namespace nb::literals;

// SeparationOracle stores unique_ptr<Separator>, so it is non-copyable.
// Tell nanobind explicitly to avoid generating copy wrappers.
namespace nanobind::detail {
template <>
struct is_copy_constructible<rcspp::sep::SeparationOracle> : std::false_type {};
}  // namespace nanobind::detail

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

    // ================================================================
    // Solver-independent algorithm API
    // ================================================================

    // --- Cut ---
    nb::class_<rcspp::sep::Cut>(m, "Cut")
        .def_prop_ro("indices", [](rcspp::sep::Cut& self) {
            return vec_view_numpy<int32_t>(self.indices, nb::find(self));
        })
        .def_prop_ro("values", [](rcspp::sep::Cut& self) {
            return vec_view_numpy<double>(self.values, nb::find(self));
        })
        .def_ro("rhs", &rcspp::sep::Cut::rhs)
        .def_ro("violation", &rcspp::sep::Cut::violation)
        .def_prop_ro("size", &rcspp::sep::Cut::size);

    // --- SeparationOracle ---
    nb::class_<rcspp::sep::SeparationOracle>(m, "SeparationOracle")
        .def(nb::init<const rcspp::Problem&>(), "problem"_a,
             nb::keep_alive<1, 2>(),  // oracle refs problem
             "Create a separation oracle for the given problem.")
        .def("add_default_separators",
             &rcspp::sep::SeparationOracle::add_default_separators,
             "Add the default separator set (SEC + RCI + Multistar + Comb).")
        .def("add_sec", [](rcspp::sep::SeparationOracle& self) {
            self.add_separator(std::make_unique<rcspp::sep::SECSeparator>());
        }, "Add Subtour Elimination Constraint separator.")
        .def("add_rci", [](rcspp::sep::SeparationOracle& self) {
            self.add_separator(std::make_unique<rcspp::sep::RCISeparator>());
        }, "Add Rounded Capacity Inequality separator.")
        .def("add_multistar", [](rcspp::sep::SeparationOracle& self) {
            self.add_separator(std::make_unique<rcspp::sep::MultistarSeparator>());
        }, "Add Multistar (GLM) separator.")
        .def("add_comb", [](rcspp::sep::SeparationOracle& self) {
            self.add_separator(std::make_unique<rcspp::sep::CombSeparator>());
        }, "Add Comb inequality separator.")
        .def("add_rglm", [](rcspp::sep::SeparationOracle& self) {
            self.add_separator(std::make_unique<rcspp::sep::RGLMSeparator>());
        }, "Add Rounded GLM separator.")
        .def("set_max_cuts_per_separator",
             &rcspp::sep::SeparationOracle::set_max_cuts_per_separator,
             "max_cuts"_a,
             "Set max cuts to keep per separator per round (0 = unlimited).")
        .def("separate", [](const rcspp::sep::SeparationOracle& self,
                            nb::ndarray<double, nb::shape<-1>> x_values,
                            nb::ndarray<double, nb::shape<-1>> y_values,
                            int32_t x_offset, int32_t y_offset,
                            double tol) {
            auto xv = x_values.view();
            auto yv = y_values.view();
            return self.separate(
                {xv.data(), static_cast<size_t>(x_values.shape(0))},
                {yv.data(), static_cast<size_t>(y_values.shape(0))},
                x_offset, y_offset, tol);
        }, "x_values"_a, "y_values"_a, "x_offset"_a = 0, "y_offset"_a = 0,
           "tol"_a = rcspp::sep::kDefaultFracTol,
           "Run separation on LP solution. Returns list of Cut objects.")
        .def("is_feasible", [](const rcspp::sep::SeparationOracle& self,
                               nb::ndarray<double, nb::shape<-1>> x_values,
                               nb::ndarray<double, nb::shape<-1>> y_values,
                               int32_t x_offset, int32_t y_offset) {
            auto xv = x_values.view();
            auto yv = y_values.view();
            return self.is_feasible(
                {xv.data(), static_cast<size_t>(x_values.shape(0))},
                {yv.data(), static_cast<size_t>(y_values.shape(0))},
                x_offset, y_offset);
        }, "x_values"_a, "y_values"_a, "x_offset"_a = 0, "y_offset"_a = 0,
           "Check if integer solution satisfies all SECs.");

    // --- WarmStartResult ---
    nb::class_<rcspp::heuristic::HeuristicResult>(m, "WarmStartResult")
        .def_prop_ro("col_values", [](rcspp::heuristic::HeuristicResult& self) {
            return vec_view_numpy<double>(self.col_values, nb::find(self));
        })
        .def_ro("objective", &rcspp::heuristic::HeuristicResult::objective);

    // --- build_warm_start ---
    m.def("build_warm_start", [](const rcspp::Problem& prob, double time_budget_ms) {
        return rcspp::heuristic::build_warm_start(prob, time_budget_ms);
    }, "problem"_a, "time_budget_ms"_a = 1000.0,
       "Build a warm-start solution via parallel randomized construction + local search.");

    // --- Preprocessing: forward/backward labeling ---
    m.def("forward_labeling", [](const rcspp::Problem& prob, int32_t source) {
        auto bounds = rcspp::preprocess::forward_labeling(prob, source);
        auto* data = new std::vector<double>(std::move(bounds));
        nb::capsule owner(data, [](void* p) noexcept {
            delete static_cast<std::vector<double>*>(p);
        });
        return nb::ndarray<nb::numpy, double, nb::shape<-1>>(
            data->data(), {data->size()}, std::move(owner));
    }, "problem"_a, "source"_a,
       "Capacity-aware 2-cycle elimination labeling from source.");

    m.def("backward_labeling", [](const rcspp::Problem& prob, int32_t target) {
        auto bounds = rcspp::preprocess::backward_labeling(prob, target);
        auto* data = new std::vector<double>(std::move(bounds));
        nb::capsule owner(data, [](void* p) noexcept {
            delete static_cast<std::vector<double>*>(p);
        });
        return nb::ndarray<nb::numpy, double, nb::shape<-1>>(
            data->data(), {data->size()}, std::move(owner));
    }, "problem"_a, "target"_a,
       "Capacity-aware 2-cycle elimination labeling from target.");

    // --- Preprocessing: edge elimination ---
    m.def("edge_elimination", [](const rcspp::Problem& prob,
                                  nb::ndarray<double, nb::shape<-1>> fwd,
                                  nb::ndarray<double, nb::shape<-1>> bwd,
                                  double upper_bound, double correction) {
        auto fv = fwd.view();
        auto bv = bwd.view();
        std::vector<double> f(fv.data(), fv.data() + fwd.shape(0));
        std::vector<double> b(bv.data(), bv.data() + bwd.shape(0));
        auto eliminated = rcspp::preprocess::edge_elimination(
            prob, f, b, upper_bound, correction);
        // Convert vector<bool> to numpy int8 array
        size_t m = eliminated.size();
        auto* data = new std::vector<int8_t>(m);
        for (size_t i = 0; i < m; ++i) (*data)[i] = eliminated[i] ? 1 : 0;
        nb::capsule owner(data, [](void* p) noexcept {
            delete static_cast<std::vector<int8_t>*>(p);
        });
        return nb::ndarray<nb::numpy, int8_t, nb::shape<-1>>(
            data->data(), {m}, std::move(owner));
    }, "problem"_a, "fwd_bounds"_a, "bwd_bounds"_a,
       "upper_bound"_a, "correction"_a,
       "Edge elimination using labeling bounds. Returns bool array per edge.");

    // --- BoundPropagator ---
    nb::class_<rcspp::preprocess::BoundPropagator>(m, "BoundPropagator")
        .def("__init__", [](rcspp::preprocess::BoundPropagator* self,
                            const rcspp::Problem& prob,
                            nb::ndarray<double, nb::shape<-1>> fwd,
                            nb::ndarray<double, nb::shape<-1>> bwd,
                            double correction) {
            auto fv = fwd.view();
            auto bv = bwd.view();
            new (self) rcspp::preprocess::BoundPropagator(
                prob,
                std::vector<double>(fv.data(), fv.data() + fwd.shape(0)),
                std::vector<double>(bv.data(), bv.data() + bwd.shape(0)),
                correction);
        }, "problem"_a, "fwd_bounds"_a, "bwd_bounds"_a, "correction"_a,
           nb::keep_alive<1, 2>())
        .def("set_all_pairs_bounds", [](rcspp::preprocess::BoundPropagator& self,
                                        nb::ndarray<double, nb::shape<-1>> dist) {
            auto dv = dist.view();
            self.set_all_pairs_bounds(
                std::vector<double>(dv.data(), dv.data() + dist.shape(0)));
        }, "dist"_a, "Set all-pairs bounds (flat n*n array) for stronger Trigger B.")
        .def("sweep", [](const rcspp::preprocess::BoundPropagator& self,
                          double upper_bound,
                          nb::ndarray<double, nb::shape<-1>> col_upper) {
            auto cv = col_upper.view();
            auto fixings = self.sweep(upper_bound,
                {cv.data(), static_cast<size_t>(col_upper.shape(0))});
            auto* data = new std::vector<int32_t>(std::move(fixings));
            nb::capsule owner(data, [](void* p) noexcept {
                delete static_cast<std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t, nb::shape<-1>>(
                data->data(), {data->size()}, std::move(owner));
        }, "upper_bound"_a, "col_upper"_a,
           "Trigger A sweep: return edge indices that can be fixed to 0.")
        .def("propagate_fixed_edge",
             [](const rcspp::preprocess::BoundPropagator& self,
                int32_t edge, double upper_bound,
                nb::ndarray<double, nb::shape<-1>> col_upper) {
            auto cv = col_upper.view();
            auto fixings = self.propagate_fixed_edge(edge, upper_bound,
                {cv.data(), static_cast<size_t>(col_upper.shape(0))});
            auto* data = new std::vector<int32_t>(std::move(fixings));
            nb::capsule owner(data, [](void* p) noexcept {
                delete static_cast<std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t, nb::shape<-1>>(
                data->data(), {data->size()}, std::move(owner));
        }, "edge"_a, "upper_bound"_a, "col_upper"_a,
           "Trigger B: given edge fixed to 1, return edges that can be fixed to 0.");

    // ================================================================
    // HiGHS integration (only when built with RCSPP_BUILD_HIGHS)
    // ================================================================
#ifdef RCSPP_HAS_HIGHS
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

    m.attr("has_highs") = true;
#else
    m.attr("has_highs") = false;
#endif

    // IO functions (always available)
    m.def("load", [](const std::string& path) {
        return rcspp::io::load(path);
    }, "path"_a, "Load an RCSPP instance from file (auto-detect format)");
}
