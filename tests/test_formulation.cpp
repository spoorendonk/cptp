#include <cstdio>
#include <filesystem>
#include <iostream>
#include <vector>

#include "core/io.h"
#include "core/problem.h"
#include "Highs.h"

// Build and solve formulation WITHOUT separators.
// Then verify a known PathWyse tour is feasible.
int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <instance.sppcc>\n", argv[0]);
        return 1;
    }

    auto prob = cptp::io::load(argv[1]);
    const auto& graph = prob.graph();
    const int32_t n = prob.num_nodes();
    const int32_t m = prob.num_edges();
    const int total_vars = m + n;

    printf("Instance: %s (%d nodes, %d edges)\n", prob.name.c_str(), n, m);
    printf("Depot: %d, Capacity: %.0f\n", prob.depot(), prob.capacity());

    // Print first few profits and costs
    printf("\nProfits (first 10):\n");
    for (int i = 0; i < std::min(n, 10); ++i) {
        printf("  node %d: profit=%.0f, demand=%.0f\n", i,
               prob.profit(i), prob.demand(i));
    }

    // Build formulation in HiGHS
    Highs highs;
    highs.setOptionValue("output_flag", false);

    // Objective: min sum(cost_e * x_e) - sum(profit_i * y_i)
    for (auto e : graph.edges()) {
        highs.addCol(prob.edge_cost(e), 0.0, 1.0, 0, nullptr, nullptr);
    }
    for (int32_t i = 0; i < n; ++i) {
        highs.addCol(-prob.profit(i), 0.0, 1.0, 0, nullptr, nullptr);
    }
    for (int i = 0; i < total_vars; ++i) {
        highs.changeColIntegrality(i, HighsVarType::kInteger);
    }
    highs.changeObjectiveSense(ObjSense::kMinimize);

    // Degree: sum_{e incident to i} x_e = 2*y_i
    for (int32_t i = 0; i < n; ++i) {
        std::vector<int> indices;
        std::vector<double> values;
        for (auto e : graph.incident_edges(i)) {
            indices.push_back(static_cast<int>(e));
            values.push_back(1.0);
        }
        indices.push_back(m + i);
        values.push_back(-2.0);
        highs.addRow(0.0, 0.0, (int)indices.size(), indices.data(), values.data());
    }

    // Depot: y_depot = 1
    {
        int idx = m + prob.depot();
        double val = 1.0;
        highs.addRow(1.0, 1.0, 1, &idx, &val);
    }

    // Capacity: sum(demand_i * y_i) <= Q
    if (prob.capacity() > 0 && prob.capacity() < 1e17) {
        std::vector<int> indices;
        std::vector<double> values;
        for (int32_t i = 0; i < n; ++i) {
            double d = prob.demand(i);
            if (d > 0) {
                indices.push_back(m + i);
                values.push_back(d);
            }
        }
        if (!indices.empty()) {
            highs.addRow(-kHighsInf, prob.capacity(),
                         (int)indices.size(), indices.data(), values.data());
        }
    }

    printf("\nSolving MIP (without SECs)...\n");
    highs.run();

    double obj_val;
    highs.getInfoValue("objective_function_value", obj_val);
    printf("HiGHS objective: %.2f (our_obj = %.2f)\n", obj_val, -obj_val);
    printf("Model status: %d\n", (int)highs.getModelStatus());

    // Print the solution tour
    const auto& sol = highs.getSolution();
    printf("\nVisited nodes (y_i > 0.5):");
    double total_demand = 0;
    for (int32_t i = 0; i < n; ++i) {
        if (sol.col_value[m + i] > 0.5) {
            printf(" %d", i);
            total_demand += prob.demand(i);
        }
    }
    printf("\nTotal demand: %.0f / %.0f\n", total_demand, prob.capacity());

    printf("\nActive edges (x_e > 0.5):\n");
    for (auto e : graph.edges()) {
        if (sol.col_value[e] > 0.5) {
            printf("  %d -- %d (cost=%.0f)\n",
                   graph.edge_source(e),
                   graph.edge_target(e),
                   prob.edge_cost(e));
        }
    }

    return 0;
}
