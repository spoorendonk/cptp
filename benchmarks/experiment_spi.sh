#!/bin/bash
# SPI (Shortest Path Inequality) separator parameter sweep
# Compares: baseline (no SPI) vs all_pairs_propagation=true (SPI + lifting)
# Across a representative set of SPPRCLIB instances.

set -euo pipefail

SOLVER="./build/cptp-solve"
INST_DIR="benchmarks/instances/spprclib"
TIME_LIMIT="${1:-60}"

# Representative mix: small (B-45..B-57), medium (A-63..B-78), hard (E-101..M-151)
INSTANCES=(
    "B-n45-k6-54.sppcc"
    "B-n50-k8-40.sppcc"
    "B-n52-k7-15.sppcc"
    "B-n57-k7-20.sppcc"
    "A-n54-k7-149.sppcc"
    "A-n60-k9-57.sppcc"
    "A-n63-k9-157.sppcc"
    "A-n64-k9-45.sppcc"
    "A-n65-k9-10.sppcc"
    "A-n69-k9-42.sppcc"
    "A-n80-k10-14.sppcc"
    "B-n66-k9-50.sppcc"
    "B-n67-k10-26.sppcc"
    "B-n68-k9-65.sppcc"
    "B-n78-k10-70.sppcc"
    "E-n76-k7-44.sppcc"
    "E-n101-k8-291.sppcc"
    "E-n101-k14-158.sppcc"
    "M-n121-k7-260.sppcc"
    "M-n151-k12-15.sppcc"
)

BASE_OPTS="--time_limit $TIME_LIMIT --output_flag false --presolve off"

parse_result() {
    local line="$1"
    echo "$line" | awk '{
        for(i=1;i<=NF;i++) {
            if($i=="Objective:") obj=$(i+1)
            if($i=="Bound:") bound=$(i+1)
            if($i=="Gap:") gap=$(i+1)
            if($i=="Time:") time=$(i+1)
            if($i=="Nodes:") nodes=$(i+1)
        }
        gsub(/%/,"",gap); gsub(/s$/,"",time)
        printf "%12s %12s %8s %8s %8s", obj, bound, gap, time, nodes
    }'
}

echo "=== SPI Separator Sweep ==="
echo "Time limit: ${TIME_LIMIT}s per run"
echo ""

printf "%-25s %4s %12s %12s %8s %8s %8s\n" "Instance" "SPI" "Objective" "Bound" "Gap%" "Time" "Nodes"
echo "---------------------------------------------------------------------------------------------"

wins_on=0
wins_off=0
ties=0

for inst in "${INSTANCES[@]}"; do
    if [[ ! -f "$INST_DIR/$inst" ]]; then
        printf "%-25s  SKIP (not found)\n" "$inst"
        continue
    fi

    # Baseline: no SPI
    result_off=$("$SOLVER" "$INST_DIR/$inst" $BASE_OPTS \
        --all_pairs_propagation false 2>/dev/null | grep "Objective:" || echo "FAILED")

    printf "%-25s %4s " "$inst" "OFF"
    if [[ "$result_off" == "FAILED" ]]; then
        echo "FAILED"
        obj_off=""
        time_off=""
    else
        parse_result "$result_off"
        echo ""
        obj_off=$(echo "$result_off" | awk '{for(i=1;i<=NF;i++) if($i=="Objective:") print $(i+1)}')
        time_off=$(echo "$result_off" | awk '{for(i=1;i<=NF;i++) if($i=="Time:") {gsub(/s$/,"",$(i+1)); print $(i+1)}}')
        gap_off=$(echo "$result_off" | awk '{for(i=1;i<=NF;i++) if($i=="Gap:") {gsub(/%$/,"",$(i+1)); print $(i+1)}}')
    fi

    # SPI enabled
    result_on=$("$SOLVER" "$INST_DIR/$inst" $BASE_OPTS \
        --all_pairs_propagation true 2>/dev/null | grep "Objective:" || echo "FAILED")

    printf "%-25s %4s " "$inst" "ON"
    if [[ "$result_on" == "FAILED" ]]; then
        echo "FAILED"
        obj_on=""
        time_on=""
    else
        parse_result "$result_on"
        echo ""
        obj_on=$(echo "$result_on" | awk '{for(i=1;i<=NF;i++) if($i=="Objective:") print $(i+1)}')
        time_on=$(echo "$result_on" | awk '{for(i=1;i<=NF;i++) if($i=="Time:") {gsub(/s$/,"",$(i+1)); print $(i+1)}}')
        gap_on=$(echo "$result_on" | awk '{for(i=1;i<=NF;i++) if($i=="Gap:") {gsub(/%$/,"",$(i+1)); print $(i+1)}}')
    fi

    # Compare
    if [[ -n "${gap_off:-}" && -n "${gap_on:-}" ]]; then
        better=$(awk "BEGIN { if ($gap_on + 0.01 < $gap_off) print \"ON\"; else if ($gap_off + 0.01 < $gap_on) print \"OFF\"; else print \"TIE\" }")
        if [[ "$better" == "ON" ]]; then
            ((wins_on++)) || true
        elif [[ "$better" == "OFF" ]]; then
            ((wins_off++)) || true
        else
            ((ties++)) || true
        fi
    fi

    echo ""
done

echo "---------------------------------------------------------------------------------------------"
echo "Summary: SPI ON wins=$wins_on | SPI OFF wins=$wins_off | Ties=$ties"
