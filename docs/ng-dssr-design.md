# ng/DSSR Design Notes

This note records the Step 0 design choices used for ng-path bounds and async
integration.

## Source inspiration

- `../blades`: bucket-graph implementations and practical dominance engineering
  (data layout and filtering strategy)
- `../pathwyse`: DSSR-style iterative ng growth driven by cycle detection

## Chosen architecture

We keep the existing branch-and-cut formulation and replace the bound engine:

1. run ng labeling with small neighborhoods
2. detect cycles in reconstructed best paths
3. grow ng sets only for nodes involved in cycles
4. rerun until convergence or `ng_max_size`

This preserves the current solver architecture while tightening bounds.

## Dominance/data-layout choices

`src/preprocess/ng_labeling.h` uses a node-local SoA-like layout:

- `cost[]`, `demand[]`, `prev_node[]`, `prev_slot[]`, `active[]`
- packed `visited[]` bitmaps for ng memory

Candidate scans are split into two passes:

- `<=` filter for "does an existing label dominate the candidate?"
- `>=` filter for "does candidate dominate existing labels?"

Both passes use AVX2 prefilters when enabled (`ng_simd=true`) and scalar
fallback otherwise.

## Async integration invariants

`src/preprocess/shared_bounds.h` defines a versioned shared snapshot.

- published bounds are monotone-tightened element-wise (`max` on finite values)
- each change increments `version`
- consumers can cheaply detect new snapshots via version change

`HiGHSBridge::install_propagator()` refreshes to the newest snapshot and applies
existing Trigger A/B logic unchanged.

## Non-goals for this slice

- no bucket-graph data structure yet
- no jump-arc machinery
- no branch-cut-price integration

The current scope focuses on tighter preprocessing/propagation bounds with low
risk to the rest of the solver.
