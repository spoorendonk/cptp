#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "sep/cut.h"
#include "sep/separation_context.h"
#include "sep/separator.h"

namespace cptp {
class Problem;
}

namespace cptp::sep {

/// Solver-independent separation oracle.
///
/// Bundles the support-graph construction, Gomory-Hu tree computation,
/// and parallel separator execution into a single call.  Designed to be
/// called from any MIP solver's cut callback.
///
/// Usage (C++):
///   SeparationOracle oracle(prob);
///   oracle.add_default_separators();
///   // inside solver callback:
///   auto cuts = oracle.separate(x_vals, y_vals, x_off, y_off);
///
/// Usage (Python — via bindings):
///   oracle = cptp_bac.SeparationOracle(prob)
///   oracle.add_default_separators()
///   cuts = oracle.separate(x_vals, y_vals, x_off, y_off)
class SeparationOracle {
 public:
    explicit SeparationOracle(const Problem& prob);

    /// Add a user-constructed separator.
    void add_separator(std::unique_ptr<Separator> sep);

    /// Add the default separator set.
    ///
    /// Tour mode (source == target): SEC + RCI + Multistar + Comb.
    /// Path mode (source != target): SEC only.
    void add_default_separators();

    /// Run all separators on an LP relaxation solution.
    ///
    /// Builds the fractional support graph, computes the Gomory-Hu tree,
    /// then executes every registered separator in parallel.
    /// Returns cuts sorted by violation (descending), limited to
    /// max_cuts_per_separator per separator.
    ///
    /// @param x_values  Edge variable values from LP (size = num_edges).
    /// @param y_values  Node variable values from LP (size = num_nodes).
    /// @param x_offset  Column offset for edge variables in the LP.
    /// @param y_offset  Column offset for node variables in the LP.
    /// @param tol       Violation tolerance (default: kDefaultFracTol).
    std::vector<Cut> separate(std::span<const double> x_values,
                              std::span<const double> y_values,
                              int32_t x_offset,
                              int32_t y_offset,
                              double tol = kDefaultFracTol) const;

    /// Check whether an integer solution satisfies all SECs.
    /// Returns true if feasible (no violated lazy constraints).
    bool is_feasible(std::span<const double> x_values,
                     std::span<const double> y_values,
                     int32_t x_offset,
                     int32_t y_offset) const;

    /// Max cuts to keep per separator per round (0 = unlimited).
    void set_max_cuts_per_separator(int32_t max_cuts) { max_cuts_per_sep_ = max_cuts; }
    int32_t max_cuts_per_separator() const { return max_cuts_per_sep_; }

    /// Access the registered separators (e.g., for statistics).
    const std::vector<std::unique_ptr<Separator>>& separators() const { return separators_; }

 private:
    const Problem& prob_;
    std::vector<std::unique_ptr<Separator>> separators_;
    int32_t max_cuts_per_sep_ = 3;
};

}  // namespace cptp::sep
