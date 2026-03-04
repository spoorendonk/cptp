#!/usr/bin/env bash
# Benchmark hyperplane branching modes
# Usage: ./benchmarks/bench_hyperplane.sh [time_limit_seconds]

set -uo pipefail

SOLVER="./build/cptp-solve"
INST_DIR="benchmarks/instances/spprclib"
TIME_LIMIT="${1:-120}"

# 14 representative instances across sizes
INSTANCES=(
    # Small (< 55 nodes)
    "B-n45-k6-54.sppcc"
    "B-n50-k8-40.sppcc"
    "B-n52-k7-15.sppcc"
    "A-n54-k7-149.sppcc"
    # Medium (55-70 nodes)
    "P-n55-k7-116.sppcc"
    "A-n63-k9-157.sppcc"
    "A-n63-k10-44.sppcc"
    "A-n65-k9-10.sppcc"
    "B-n66-k9-50.sppcc"
    "B-n68-k9-65.sppcc"
    # Large (76+ nodes)
    "E-n76-k7-44.sppcc"
    "A-n80-k10-14.sppcc"
    "E-n101-k8-291.sppcc"
    "E-n101-k14-158.sppcc"
)

MODES=("off" "pairs" "clusters" "demand" "all")

# Check solver exists
if [ ! -x "$SOLVER" ]; then
    echo "ERROR: $SOLVER not found. Build first." >&2
    exit 1
fi

# Header
printf "%-28s %-8s %10s %10s %8s %8s %8s\n" \
    "Instance" "Mode" "Objective" "Bound" "Gap%" "Time(s)" "Nodes"
printf "%s\n" "$(printf '=%.0s' {1..90})"

for inst in "${INSTANCES[@]}"; do
    path="$INST_DIR/$inst"
    if [ ! -f "$path" ]; then
        echo "SKIP: $path not found" >&2
        continue
    fi
    for mode in "${MODES[@]}"; do
        # Run solver, capture the Objective summary line
        output=$("$SOLVER" "$path" \
            --time_limit "$TIME_LIMIT" \
            --output_flag false \
            --branch_hyper "$mode" 2>&1 | grep "^Objective:" || true)

        if [ -z "$output" ]; then
            printf "%-28s %-8s %10s %10s %7s%% %8s %8s\n" \
                "$inst" "$mode" "ERR" "ERR" "ERR" "ERR" "ERR"
            continue
        fi

        # Parse: "Objective: -74278  Bound: -74278  Gap: 0%  Time: 1.23s  Nodes: 45"
        obj=$(echo "$output" | grep -oP 'Objective:\s*\K[-0-9.e+]+' || echo "?")
        bound=$(echo "$output" | grep -oP 'Bound:\s*\K[-0-9.e+]+' || echo "?")
        gap=$(echo "$output" | grep -oP 'Gap:\s*\K[0-9.e+]+' || echo "?")
        time=$(echo "$output" | grep -oP 'Time:\s*\K[0-9.]+' || echo "?")
        nodes=$(echo "$output" | grep -oP 'Nodes:\s*\K[0-9]+' || echo "?")

        printf "%-28s %-8s %10s %10s %7s%% %8s %8s\n" \
            "$inst" "$mode" "$obj" "$bound" "$gap" "$time" "$nodes"
    done
done
