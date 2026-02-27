# ng/DSSR Bounds Implementation

This page documents the implemented ng-path / DSSR pipeline used in
preprocessing and bound tightening.

## What is implemented

The solver keeps the branch-and-cut core unchanged and strengthens bounds by
iterative ng growth:

1. run ng labeling with a small neighborhood size
2. reconstruct candidate paths and detect cycles
3. enlarge ng sets only at cycle-critical vertices
4. rerun labeling until convergence or `ng_max_size`

This is implemented in the preprocessing stack and consumed by model-side
propagation/fixing callbacks.

## Label storage and dominance

`src/preprocess/ng_labeling.h` uses a node-local label storage layout with
separate arrays for label attributes and compact ng-memory bitmaps. The layout
is optimized for dominance scans.

Dominance is evaluated in two phases:

- existing-label precheck: can any active label dominate the candidate?
- candidate sweep: which existing labels are dominated by the candidate?

When enabled (`ng_simd=true`), AVX2 prefilters are used before scalar checks.
The scalar path remains the reference fallback.

## Shared bounds and asynchronous handoff

`src/preprocess/shared_bounds.h` provides a versioned snapshot container for
bounds produced by ng/DSSR.

Implemented invariants:

- bounds are published with monotone tightening on finite entries
- each publish increments `version`
- consumers use version checks for cheap freshness detection

`src/model/highs_bridge.cpp` refreshes from the latest published snapshot and
applies standard Trigger A/B propagation/fixing logic without changing the LP
formulation.

## Incumbent/proof integration

The ng/DSSR pipeline can hand off:

- feasible complete paths as incumbent candidates
- proof conditions that allow early optimal termination when satisfied

This integration is active in the HiGHS bridge callback flow.

## Current scope

Implemented and enabled:

- iterative ng growth for tighter bounds
- shared snapshot publication/consumption for asynchronous updates
- reduced-cost fixing policy support (`adaptive` default)

Not implemented in this module:

- bucket-graph data structures
- jump-arc machinery
- branch-cut-price architecture changes
