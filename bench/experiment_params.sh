#!/bin/bash
# Experiment 1+2: Strong branching (mip_pscost_minreliable) and heuristic effort sweep
# Run on the 3 hardest SPPRCLIB instances with 120s time limit

set -euo pipefail

SOLVER="./build/cptp-solve"
INST_DIR="bench/instances/spprclib"
TIME_LIMIT=120

INSTANCES=(
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

echo "=== Parameter Sweep Experiment ==="
echo "Time limit: ${TIME_LIMIT}s per run"
echo ""

# Experiment 1: mip_pscost_minreliable sweep
echo "--- Experiment 1: mip_pscost_minreliable ---"
printf "%-30s %6s %12s %12s %8s %8s %8s\n" "Instance" "PsCost" "Objective" "Bound" "Gap%" "Time" "Nodes"
echo "--------------------------------------------------------------------------------------------"

for inst in "${INSTANCES[@]}"; do
    for pscost in 0 2 4 8 16; do
        result_line=$("$SOLVER" "$INST_DIR/$inst" $BASE_OPTS \
            --mip_pscost_minreliable "$pscost" 2>/dev/null | grep "Objective:" || echo "FAILED")

        printf "%-30s %6d " "$inst" "$pscost"
        if [[ "$result_line" == "FAILED" ]]; then
            echo "FAILED"
        else
            parse_result "$result_line"
            echo ""
        fi
    done
    echo ""
done

# Experiment 2: mip_heuristic_effort sweep (with default pscost=8)
echo "--- Experiment 2: mip_heuristic_effort ---"
printf "%-30s %8s %12s %12s %8s %8s %8s\n" "Instance" "HeurEff" "Objective" "Bound" "Gap%" "Time" "Nodes"
echo "--------------------------------------------------------------------------------------------"

for inst in "${INSTANCES[@]}"; do
    for heur_eff in 0 0.01 0.05 0.1; do
        result_line=$("$SOLVER" "$INST_DIR/$inst" $BASE_OPTS \
            --mip_heuristic_effort "$heur_eff" 2>/dev/null | grep "Objective:" || echo "FAILED")

        printf "%-30s %8s " "$inst" "$heur_eff"
        if [[ "$result_line" == "FAILED" ]]; then
            echo "FAILED"
        else
            parse_result "$result_line"
            echo ""
        fi
    done
    echo ""
done
