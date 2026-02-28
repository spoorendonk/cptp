# Solve Concurrency

File: `src/model/model.cpp` (`SolveConcurrencyGuard` in anonymous namespace)

## Overview

`max_concurrent_solves` controls in-process admission for concurrent
`Model::solve()` calls.

- `0` (default): this solve requests no explicit cap
- `N > 0`: this solve requests at most `N` concurrent solves

This mechanism is process-local and independent of HiGHS internal
parallelism (`--threads`).

## Effective Cap Rule

The guard tracks:

- active solve count
- multiset of finite requested caps from active and waiting solves

At admission time, a solve waits until:

- `active_solves < min(requested_limit, current_global_cap)`

where `current_global_cap` is the minimum finite cap currently registered.

Consequence: stricter requests immediately constrain later arrivals.

## Why pending limits are registered early

Finite requested caps are inserted before waiting. This prevents a race where a
strict request arrives but looser calls continue entering ahead of it.

## Flow

```mermaid
flowchart TD
    A[Construct SolveConcurrencyGuard(limit)] --> B{limit > 0}
    B -->|yes| C[Register finite cap in shared multiset]
    B -->|no| D[Use no-cap sentinel]
    C --> E[Wait until active < min(requested, current_cap)]
    D --> E
    E --> F[Increment active_solves]
    F --> G[Run solve]
    G --> H[Destructor]
    H --> I{Had finite cap}
    I -->|yes| J[Remove cap from multiset]
    I -->|no| K[Skip removal]
    J --> L[Decrement active_solves and notify all]
    K --> L
```

## Pseudocode

```text
constructor(limit):
    requested <- (limit > 0) ? limit : INF
    if requested is finite:
        insert requested into shared active_limits

    wait until active_solves < min(requested, min(active_limits or INF))
    active_solves += 1

destructor():
    if requested is finite:
        erase one requested from active_limits
    active_solves -= 1
    notify_all waiting solves
```

## Interaction with `threads`

- `max_concurrent_solves`: number of solves admitted concurrently in one process
- `threads`: worker threads used inside one admitted HiGHS solve

Tuning both is workload-dependent. Example: on shared machines, lower
`max_concurrent_solves` can reduce oversubscription even if `threads` stays high.

## Usage

```bash
./build/rcspp-solve instance.sppcc --max_concurrent_solves 1 --threads 8
```
