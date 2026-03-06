#!/usr/bin/env bash
# run_benchmarks.sh — Run cptp-solve on benchmark instances and produce CSV results.
#
# Runs build/cptp-solve on all spprclib + roberti instances, parses stdout,
# and writes rows to benchmarks/cptp.csv (one row per instance, replaced on re-run).
#
# Usage:
#   ./benchmarks/run_benchmarks.sh [--time-limit N] [instance_or_dir]
#
# Arguments:
#   --time-limit N       Per-instance time limit in seconds (default: 3600)
#   instance_or_dir      Single instance file or directory of instances.
#                        Default: both benchmarks/instances/spprclib and
#                        benchmarks/instances/roberti.
#
# Output:
#   benchmarks/cptp.csv — committed CSV with timestamp and machine columns.
#   Existing rows for re-run instances are replaced; other rows are preserved.
#
# Parsed fields from cptp-solve stdout:
#   Instance: <name> (<N> nodes, <E> edges, ...)    → instance, nodes_graph, edges
#   Objective: <obj>  Bound: <bound>  Gap: <gap>%   → obj, bound, gap_pct
#     Time: <T>s  Nodes: <N>                        → time_s, bb_nodes
#   User cuts: <N> (<M> rounds)                     → total_cuts, cut_rounds
#   <SEP>  <N> cuts  <M> rounds  <T>s               → per-separator stats
#   Local search last progress row                   → warmstart_ub, warmstart_time
#
# Prerequisites:
#   cptp-solve must be built: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SOLVER="$REPO_DIR/build/cptp-solve"

TIME_LIMIT=3600
TARGETS=()

usage() {
    sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-limit)
            TIME_LIMIT="$2"; shift 2 ;;
        --help|-h)
            usage 0 ;;
        -*)
            echo "Unknown option: $1" >&2; usage 1 ;;
        *)
            TARGETS+=("$1"); shift ;;
    esac
done

if [[ ! -x "$SOLVER" ]]; then
    echo "Error: cptp-solve not found at $SOLVER" >&2
    echo "Build first: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build" >&2
    exit 1
fi

# Collect instance files
collect_instances() {
    local target="$1"
    if [[ -f "$target" ]]; then
        echo "$target"
    elif [[ -d "$target" ]]; then
        find "$target" -maxdepth 1 \( -name '*.sppcc' -o -name '*.vrp' \) | sort
    else
        echo "Error: $target not found" >&2
        exit 1
    fi
}

INSTANCES=()
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=(
        "$SCRIPT_DIR/instances/spprclib"
        "$SCRIPT_DIR/instances/roberti"
    )
fi
for target in "${TARGETS[@]}"; do
    while IFS= read -r f; do
        INSTANCES+=("$f")
    done < <(collect_instances "$target")
done

if [[ ${#INSTANCES[@]} -eq 0 ]]; then
    echo "No .sppcc or .vrp instance files found" >&2
    exit 1
fi

OUTFILE="$SCRIPT_DIR/cptp.csv"
HEADER="instance,nodes_graph,edges,obj,bound,gap_pct,time_s,time_limit,bb_nodes,lp_iters,lp_iters_strongbr,lp_iters_sep,lp_iters_heur,total_cuts,cut_rounds,sec_cuts,sec_rounds,sec_time,rci_cuts,rci_rounds,rci_time,multistar_cuts,multistar_rounds,multistar_time,rglm_cuts,rglm_rounds,rglm_time,comb_cuts,comb_rounds,comb_time,warmstart_ub,warmstart_time,timestamp,machine"

# Create CSV with header if it doesn't exist
if [[ ! -f "$OUTFILE" ]]; then
    echo "$HEADER" > "$OUTFILE"
fi

# Detect machine info
MACHINE="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g') $(nproc)T"

# Collect instance stems to be (re-)run, then remove their existing rows
INSTANCE_STEMS=()
for inst in "${INSTANCES[@]}"; do
    INSTANCE_STEMS+=("$(basename "${inst%.*}")")
done

# Build a temporary file without rows for instances we're about to run
TMPFILE="$(mktemp)"
head -1 "$OUTFILE" > "$TMPFILE"
tail -n +2 "$OUTFILE" | while IFS= read -r line; do
    row_instance="${line%%,*}"
    skip=false
    for stem in "${INSTANCE_STEMS[@]}"; do
        if [[ "$row_instance" == "$stem" ]]; then
            skip=true
            break
        fi
    done
    if ! $skip; then
        echo "$line" >> "$TMPFILE"
    fi
done
cp "$TMPFILE" "$OUTFILE"
rm -f "$TMPFILE"

echo "cptp benchmark run — ${#INSTANCES[@]} instances, time limit: ${TIME_LIMIT}s"
echo "Output: $OUTFILE"
echo
printf "%-35s %12s %10s %8s %7s\n" "Instance" "Obj" "Time" "Nodes" "Cuts"
printf '%s\n' "$(printf '%.0s-' {1..75})"

for inst in "${INSTANCES[@]}"; do
    name="$(basename "${inst%.*}")"

    # Run solver, capture stdout
    # Keep output_flag on (default) so we can parse warmstart/preprocess logs.
    # HiGHS verbose output is interleaved but our grep patterns are specific enough.
    output="$("$SOLVER" "$inst" --time_limit "$TIME_LIMIT" 2>&1)" || true

    # Parse Instance line
    nodes_graph="" edges=""
    if inst_line="$(echo "$output" | grep -oP 'Instance: .* \(\K[0-9]+ nodes, [0-9]+ edges')"; then
        nodes_graph="$(echo "$inst_line" | grep -oP '^[0-9]+')"
        edges="$(echo "$inst_line" | grep -oP '[0-9]+(?= edges)')"
    fi

    # Parse Objective line
    obj="" bound="" gap_pct="" time_s="" bb_nodes=""
    if obj_line="$(echo "$output" | grep -P '^Objective:')"; then
        obj="$(echo "$obj_line" | grep -oP 'Objective: \K[-0-9.e+]+')"
        bound="$(echo "$obj_line" | grep -oP 'Bound: \K[-0-9.e+]+')"
        gap_pct="$(echo "$obj_line" | grep -oP 'Gap: \K[0-9.e+-]+')"
        time_s="$(echo "$obj_line" | grep -oP 'Time: \K[0-9.]+')"
        bb_nodes="$(echo "$obj_line" | grep -oP 'Nodes: \K[0-9]+')"
    fi

    # Parse User cuts
    total_cuts="" cut_rounds=""
    if cuts_line="$(echo "$output" | grep -P '^User cuts:')"; then
        total_cuts="$(echo "$cuts_line" | grep -oP 'User cuts: \K[0-9]+')"
        cut_rounds="$(echo "$cuts_line" | grep -oP '\(\K[0-9]+')"
    fi

    # Parse per-separator stats
    parse_sep() {
        local sep_name="$1"
        local line
        line="$(echo "$output" | grep -P "^\s+${sep_name}\s" | head -1)" || true
        if [[ -n "$line" ]]; then
            local cuts rounds stime
            cuts="$(echo "$line" | grep -oP '[0-9]+(?= cuts)')" || cuts=""
            rounds="$(echo "$line" | grep -oP '[0-9]+(?= rounds)')" || rounds=""
            stime="$(echo "$line" | grep -oP '[0-9.]+(?=s\s*$)')" || stime=""
            echo "${cuts},${rounds},${stime}"
        else
            echo ",,"
        fi
    }

    sec_stats="$(parse_sep SEC)"
    rci_stats="$(parse_sep RCI)"
    multistar_stats="$(parse_sep Multistar)"
    rglm_stats="$(parse_sep RGLM)"
    comb_stats="$(parse_sep Comb)"

    # Parse warmstart: extract lines between "Local search:" and "Preprocess restart:"
    # The progress rows have format: "  <starts> <iter_accum> <ub> <impr> <time>s"
    warmstart_ub="" warmstart_time=""
    ws_block="$(echo "$output" | sed -n '/^Local search:/,/^Preprocess restart:/p')" || true
    ws_line="$(echo "$ws_block" | grep -P '^\s+\d+\s+\d+\s+[-0-9.e+]+\s+\d+\s+[0-9.]+s\s*$' | tail -1)" || true
    if [[ -n "$ws_line" ]]; then
        warmstart_ub="$(echo "$ws_line" | awk '{print $3}')"
        warmstart_time="$(echo "$ws_line" | awk '{gsub(/s$/,"",$5); print $5}')"
    fi

    # lp_iters fields are not printed by cptp-solve stdout; leave blank
    lp_iters="" lp_iters_strongbr="" lp_iters_sep="" lp_iters_heur=""

    # Write CSV row with timestamp and machine
    ROW_TIMESTAMP="$(date -Iseconds)"
    echo "${name},${nodes_graph},${edges},${obj},${bound},${gap_pct},${time_s},${TIME_LIMIT},${bb_nodes},${lp_iters},${lp_iters_strongbr},${lp_iters_sep},${lp_iters_heur},${total_cuts},${cut_rounds},${sec_stats},${rci_stats},${multistar_stats},${rglm_stats},${comb_stats},${warmstart_ub},${warmstart_time},${ROW_TIMESTAMP},\"${MACHINE}\"" >> "$OUTFILE"

    # Progress line
    printf "%-35s %12s %10s %8s %7s\n" "$name" "${obj:-N/A}" "${time_s:+${time_s}s}" "${bb_nodes:-N/A}" "${total_cuts:-N/A}"
done

printf '%s\n' "$(printf '%.0s-' {1..75})"
echo "Done. Results written to $OUTFILE"
