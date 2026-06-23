#!/usr/bin/env bash
# run_ablation.sh — per-component ablation for the cptp branch-and-cut.
#
# Runs build/cptp-solve over the benchmark instances under several cut/fixing
# configurations and writes one CSV row per (config, instance) to ablation.csv.
# The contrast between configs a_sec and b_cap identifies the capacity-binding
# instances on which capacity-class cuts pay off; spi_off vs spi_on measures the
# effect of the SPI separator.
#
# Usage:
#   ./benchmarks/run_ablation.sh [--time-limit N] [--configs "a_sec b_cap ..."] [instance_or_dir ...]
#
#   --time-limit N    Per-instance wall-clock limit in seconds (default: 3600).
#   --configs "..."   Space-separated subset of config names (default: all below).
#   instance_or_dir   Instance file(s)/dir(s). Default: benchmarks/instances/{spprclib,roberti}.
#
# Output: benchmarks/ablation.csv. Existing rows for the same (config, instance)
#   pairs are replaced; other rows are preserved.
#
# Config table (name : extra flags passed to cptp-solve):
#   a_sec        SEC only                         --enable_rci false --enable_multistar false
#   b_cap        + capacity cuts (default set)    (no extra flags: SEC+RCI+Multistar)
#   c_combrglm   + comb + rounded GLM             --enable_comb true --enable_rglm true
#   d_nofix      c, fixing/propagation OFF        --rc_fixing off --edge_elimination false
#                                                 --bounds_propagation false (+ comb/rglm)
#   d_fix        c, fixing/propagation ON         --rc_fixing adaptive (+ comb/rglm)
#   spi_off      default set + all-pairs, no SPI  --all_pairs_bounds true
#   spi_on       default set + all-pairs + SPI    --all_pairs_bounds true --enable_spi true
#
# Prerequisite: cptp-solve must be built (cmake -B build && cmake --build build).

set -euo pipefail

round3() { [[ -n "${1:-}" ]] && awk -v val="$1" 'BEGIN{v=val+0; if(v>999999||v<-999999){print val}else{x=sprintf("%.3f",v); sub(/\.?0+$/,"",x); print x}}' || echo ""; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SOLVER="$REPO_DIR/build/cptp-solve"
TIME_LIMIT=3600
CONFIGS_REQUESTED=""
TARGETS=()

config_flags() {
    case "$1" in
        a_sec)      echo "--enable_rci false --enable_multistar false" ;;
        b_cap)      echo "" ;;
        c_combrglm) echo "--enable_comb true --enable_rglm true" ;;
        d_nofix)    echo "--enable_comb true --enable_rglm true --rc_fixing off --edge_elimination false --bounds_propagation false" ;;
        d_fix)      echo "--enable_comb true --enable_rglm true --rc_fixing adaptive" ;;
        spi_off)    echo "--all_pairs_bounds true" ;;
        spi_on)     echo "--all_pairs_bounds true --enable_spi true" ;;
        *)          echo "__UNKNOWN__" ;;
    esac
}
ALL_CONFIGS="a_sec b_cap c_combrglm d_nofix d_fix spi_off spi_on"

usage() { sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-limit) TIME_LIMIT="$2"; shift 2 ;;
        --configs)    CONFIGS_REQUESTED="$2"; shift 2 ;;
        --help|-h)    usage 0 ;;
        -*)           echo "Unknown option: $1" >&2; usage 1 ;;
        *)            TARGETS+=("$1"); shift ;;
    esac
done

if [[ ! -x "$SOLVER" ]]; then
    echo "Error: cptp-solve not found at $SOLVER" >&2
    echo "Build first: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build" >&2
    exit 1
fi

CONFIGS="${CONFIGS_REQUESTED:-$ALL_CONFIGS}"
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=("$SCRIPT_DIR/instances/spprclib" "$SCRIPT_DIR/instances/roberti")
fi

collect_instances() {
    local target="$1"
    if [[ -f "$target" ]]; then echo "$target"
    elif [[ -d "$target" ]]; then find "$target" -maxdepth 1 \( -name '*.sppcc' -o -name '*.vrp' \) | sort
    else echo "Error: $target not found" >&2; exit 1; fi
}

INSTANCES=()
for target in "${TARGETS[@]}"; do
    while IFS= read -r f; do INSTANCES+=("$f"); done < <(collect_instances "$target")
done
[[ ${#INSTANCES[@]} -gt 0 ]] || { echo "No .sppcc or .vrp instance files found" >&2; exit 1; }

OUTFILE="$SCRIPT_DIR/ablation.csv"
HEADER="config,instance,nodes_graph,edges,obj,bound,gap_pct,time_s,time_limit,bb_nodes,total_cuts,cut_rounds,sec_cuts,sec_rounds,sec_time,rci_cuts,rci_rounds,rci_time,multistar_cuts,multistar_rounds,multistar_time,rglm_cuts,rglm_rounds,rglm_time,comb_cuts,comb_rounds,comb_time,spi_cuts,spi_rounds,spi_time,timestamp,machine"
[[ -f "$OUTFILE" ]] || echo "$HEADER" > "$OUTFILE"

MACHINE="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g') $(nproc)T"

# Drop existing rows for the (config,instance) pairs we are about to (re-)run.
declare -A RERUN
for cfg in $CONFIGS; do
    for inst in "${INSTANCES[@]}"; do RERUN["$cfg,$(basename "${inst%.*}")"]=1; done
done
TMPFILE="$(mktemp)"; trap 'rm -f "$TMPFILE"' EXIT
head -1 "$OUTFILE" > "$TMPFILE"
tail -n +2 "$OUTFILE" | while IFS= read -r line; do
    key="$(echo "$line" | cut -d, -f1,2)"
    [[ -n "${RERUN[$key]:-}" ]] || echo "$line" >> "$TMPFILE"
done
cp "$TMPFILE" "$OUTFILE"; rm -f "$TMPFILE"; trap - EXIT

echo "cptp ablation — ${#INSTANCES[@]} instances × configs: $CONFIGS — time limit ${TIME_LIMIT}s"
echo "Output: $OUTFILE"
printf "%-10s %-30s %12s %9s %8s\n" "Config" "Instance" "Obj" "Time" "Nodes"
printf '%s\n' "$(printf '%.0s-' {1..72})"

parse_sep() {
    local line; line="$(echo "$OUTPUT" | grep -P "^\s+$1\s" | head -1)" || true
    if [[ -n "$line" ]]; then
        local c r t
        c="$(echo "$line" | grep -oP '[0-9]+(?= cuts)')" || c=""
        r="$(echo "$line" | grep -oP '[0-9]+(?= rounds)')" || r=""
        t="$(round3 "$(echo "$line" | grep -oP '[0-9.]+(?=s\s*$)')")" || t=""
        echo "${c},${r},${t}"
    else echo ",,"; fi
}

for cfg in $CONFIGS; do
    flags="$(config_flags "$cfg")"
    [[ "$flags" != "__UNKNOWN__" ]] || { echo "Unknown config: $cfg" >&2; exit 1; }
    for inst in "${INSTANCES[@]}"; do
        name="$(basename "${inst%.*}")"
        # shellcheck disable=SC2086
        OUTPUT="$("$SOLVER" "$inst" --time_limit "$TIME_LIMIT" $flags 2>&1)" || true

        nodes_graph="" edges=""
        if il="$(echo "$OUTPUT" | grep -oP 'Instance: .* \(\K[0-9]+ nodes, [0-9]+ edges')"; then
            nodes_graph="$(echo "$il" | grep -oP '^[0-9]+')"
            edges="$(echo "$il" | grep -oP '[0-9]+(?= edges)')"
        fi
        obj="" bound="" gap_pct="" time_s="" bb_nodes=""
        if ol="$(echo "$OUTPUT" | grep -P '^Objective:')"; then
            obj="$(round3 "$(echo "$ol" | grep -oP 'Objective: \K[-0-9.e+]+')")"
            bound="$(round3 "$(echo "$ol" | grep -oP 'Bound: \K[-0-9.e+]+')")"
            gap_pct="$(round3 "$(echo "$ol" | grep -oP 'Gap: \K[0-9.e+-]+')")"
            time_s="$(round3 "$(echo "$ol" | grep -oP 'Time: \K[0-9.]+')")"
            bb_nodes="$(echo "$ol" | grep -oP 'Nodes: \K[0-9]+')"
        fi
        total_cuts="" cut_rounds=""
        if cl="$(echo "$OUTPUT" | grep -P '^User cuts:')"; then
            total_cuts="$(echo "$cl" | grep -oP 'User cuts: \K[0-9]+')"
            cut_rounds="$(echo "$cl" | grep -oP '\(\K[0-9]+')"
        fi
        sec_stats="$(parse_sep SEC)"; rci_stats="$(parse_sep RCI)"
        multistar_stats="$(parse_sep Multistar)"; rglm_stats="$(parse_sep RGLM)"
        comb_stats="$(parse_sep Comb)"; spi_stats="$(parse_sep SPI)"
        ts="$(date -Iseconds)"

        echo "${cfg},${name},${nodes_graph},${edges},${obj},${bound},${gap_pct},${time_s},${TIME_LIMIT},${bb_nodes},${total_cuts},${cut_rounds},${sec_stats},${rci_stats},${multistar_stats},${rglm_stats},${comb_stats},${spi_stats},${ts},\"${MACHINE}\"" >> "$OUTFILE"
        printf "%-10s %-30s %12s %9s %8s\n" "$cfg" "$name" "${obj:-NA}" "${time_s:-NA}" "${bb_nodes:-NA}"
    done
done

echo "Done. Rows written to $OUTFILE"
