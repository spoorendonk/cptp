#!/usr/bin/env bash
# run_benchmarks.sh — Run cptp-solve on benchmark instances and produce CSV results.
#
# Runs build/cptp-solve on all spprclib + roberti instances, parses stdout,
# and writes rows to benchmarks/cptp.csv (one row per instance, replaced on re-run).
#
# Usage:
#   ./benchmarks/run_benchmarks.sh [--time-limit N] [--max-hours H] [--resume] [instance_or_dir]
#
# Arguments:
#   --time-limit N       Per-instance time limit in seconds (default: 3600)
#   --max-hours H        Overall wall-clock budget for THIS invocation in hours
#                        (float, default: 0 = unlimited). Run the sweep "H hours
#                        at a time": stops once H hours elapse and skips already-
#                        completed rows, so re-running the same command continues
#                        where it left off. Implies --resume.
#   --resume             Skip instances already present in the CSV and run only
#                        the missing ones (no clock). Without --max-hours/--resume,
#                        existing rows for the requested instances are replaced.
#   instance_or_dir      Single instance file or directory of instances.
#                        Default: both benchmarks/instances/spprclib and
#                        benchmarks/instances/roberti.
#
# Output:
#   benchmarks/cptp.csv — committed CSV with timestamp and machine columns.
#   benchmarks/logs/<suite>/<instance>.log.gz — full per-run stdout (gitignored).
#   Default mode replaces rows for the requested instances; --resume/--max-hours
#   skip already-completed rows instead. If the CSV header is out of date (e.g.
#   after a schema change), the old file is archived to cptp.csv.bak first.
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

# Round a number to at most 3 decimal places, stripping trailing zeros.
round3() { [[ -n "${1:-}" ]] && awk -v val="$1" 'BEGIN{v=val+0; if(v>999999||v<-999999){print val}else{x=sprintf("%.3f",v); sub(/\.?0+$/,"",x); print x}}' || echo ""; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SOLVER="$REPO_DIR/build/cptp-solve"

TIME_LIMIT=3600
MAX_HOURS=0
RESUME=0
TARGETS=()

usage() {
    sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-limit)
            TIME_LIMIT="$2"; shift 2 ;;
        --max-hours)
            MAX_HOURS="$2"; shift 2 ;;
        --resume)
            RESUME=1; shift ;;
        --help|-h)
            usage 0 ;;
        -*)
            echo "Unknown option: $1" >&2; usage 1 ;;
        *)
            TARGETS+=("$1"); shift ;;
    esac
done

# Wall-clock budget (seconds) for this invocation; 0 = unlimited. A budget
# implies resume/skip-completed so repeated "H hours at a time" runs progress.
# Validate up front: a non-numeric/negative value must error, not silently
# parse to 0 (= unlimited) and run unbounded.
if ! [[ "$MAX_HOURS" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: --max-hours must be a non-negative number (got '$MAX_HOURS')" >&2
    exit 1
fi
BUDGET_SEC="$(awk -v h="$MAX_HOURS" 'BEGIN{printf "%d", h*3600}')"
RESUME_SKIP=0
if [[ "$RESUME" -eq 1 || "$BUDGET_SEC" -gt 0 ]]; then RESUME_SKIP=1; fi

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
LOGROOT="$SCRIPT_DIR/logs"
# New phase/fixing columns are appended after `machine` (append-only: every
# existing column keeps its index; paper scripts key by name). sep_time is the
# separation-phase wall-clock (one span per round around the parallel separators,
# so it credits the parallelism); prop_time includes rcfix_time (RC fixing runs
# inside the propagator); heur_time is the LP-guided B&C callback only (the
# warm-start heuristic is the separate warmstart_time column).
HEADER="instance,nodes_graph,edges,obj,bound,gap_pct,time_s,time_limit,bb_nodes,lp_iters,lp_iters_strongbr,lp_iters_sep,lp_iters_heur,total_cuts,cut_rounds,sec_cuts,sec_rounds,sec_time,rci_cuts,rci_rounds,rci_time,multistar_cuts,multistar_rounds,multistar_time,rglm_cuts,rglm_rounds,rglm_time,comb_cuts,comb_rounds,comb_time,warmstart_ub,warmstart_time,timestamp,machine,sep_time,prop_time,heur_time,rcfix_time,sweep_fixings,chain_fixings,node_fixings,rc_fix0,rc_fix1"

# Create the CSV, or migrate a stale-schema file. The schema change here is
# append-only, so widen existing rows in place — pad each with the new trailing
# empty fields — instead of archiving. This preserves all committed data and
# keeps subset/--resume runs safe (widened rows have empty new columns, so the
# done-set below treats them as not-yet-run). A non-append-only header change
# (reordered/removed columns) falls back to a timestamped archive.
if [[ ! -f "$OUTFILE" ]]; then
    echo "$HEADER" > "$OUTFILE"
else
    OLD_HEADER="$(head -1 "$OUTFILE")"
    if [[ "$OLD_HEADER" != "$HEADER" ]]; then
        if [[ "$HEADER" == "$OLD_HEADER",* ]]; then
            old_n="$(awk -F, '{print NF; exit}' <<<"$OLD_HEADER")"
            new_n="$(awk -F, '{print NF; exit}' <<<"$HEADER")"
            pad="$(printf ',%.0s' $(seq 1 $((new_n - old_n))))"
            TMPMIG="$(mktemp)"
            { echo "$HEADER"; tail -n +2 "$OUTFILE" | sed "s/\$/${pad}/"; } > "$TMPMIG"
            mv "$TMPMIG" "$OUTFILE"
            echo "Widened $OUTFILE to new schema (+$((new_n - old_n)) columns; existing rows padded)." >&2
        else
            BAK="$OUTFILE.$(date +%Y%m%d%H%M%S).bak"
            echo "CSV header incompatible with new schema; archiving old file -> $BAK" >&2
            mv "$OUTFILE" "$BAK"
            echo "$HEADER" > "$OUTFILE"
        fi
    fi
fi

# Detect machine info
MACHINE="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g') $(nproc)T"

# Collect instance stems to be (re-)run.
INSTANCE_STEMS=()
for inst in "${INSTANCES[@]}"; do
    INSTANCE_STEMS+=("$(basename "${inst%.*}")")
done

# Drop any existing CSV row for an instance (exact first-field match), so a fresh
# run replaces it rather than appending a duplicate. Called per-instance right
# before a run, so a budget-stop leaves not-yet-run rows (incl. widened/padded
# ones) intact. Atomic via temp + mv.
drop_row() {
    local key="$1" tmp
    tmp="$(mktemp)"
    awk -F, -v k="$key" 'NR==1 || $1!=k' "$OUTFILE" > "$tmp"
    mv "$tmp" "$OUTFILE"
}

# Completed-instance set (for resume/skip).
# A row counts as completed only if its new columns are populated (last field,
# rc_fix1, non-empty); rows widened from an old schema have empty new fields, so
# resume re-runs them to fill the new columns instead of skipping them.
declare -A DONE
while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    [[ -n "${line##*,}" ]] || continue
    DONE["${line%%,*}"]=1
done < <(tail -n +2 "$OUTFILE")

RUN_START="$(date +%s)"
ran_count=0
skipped_count=0
stopped_early=0

echo "cptp benchmark run — ${#INSTANCES[@]} instances, time limit: ${TIME_LIMIT}s"
echo "Output: $OUTFILE"
echo
printf "%-35s %12s %10s %8s %7s\n" "Instance" "Obj" "Time" "Nodes" "Cuts"
printf '%s\n' "$(printf '%.0s-' {1..75})"

for inst in "${INSTANCES[@]}"; do
    name="$(basename "${inst%.*}")"

    # Stop before starting a new instance once the wall budget is spent. An
    # already-running instance still finishes (overshoot <= one --time-limit).
    if [[ "$BUDGET_SEC" -gt 0 && $(( $(date +%s) - RUN_START )) -ge "$BUDGET_SEC" ]]; then
        echo "Wall budget ${MAX_HOURS}h reached — stopping; re-run to continue."
        stopped_early=1
        break
    fi

    # Resume/skip: leave already-completed rows untouched.
    if [[ "$RESUME_SKIP" -eq 1 && -n "${DONE[$name]:-}" ]]; then
        skipped_count=$((skipped_count + 1))
        continue
    fi

    # Replace any existing row for this instance (stale, or widened/padded).
    drop_row "$name"

    # Run solver; persist full stdout (suite = instance's parent dir), parse from
    # the saved file so the CSV and log always agree, then gzip it.
    # Keep output_flag on (default) so we can parse warmstart/preprocess logs.
    suite="$(basename "$(dirname "$inst")")"
    logdir="$LOGROOT/$suite"
    mkdir -p "$logdir"
    logfile="$logdir/$name.log"
    "$SOLVER" "$inst" --time_limit "$TIME_LIMIT" > "$logfile" 2>&1 || true
    output="$(cat "$logfile")"
    gzip -f "$logfile"

    # Parse Instance line
    nodes_graph="" edges=""
    if inst_line="$(echo "$output" | grep -oP 'Instance: .* \(\K[0-9]+ nodes, [0-9]+ edges')"; then
        nodes_graph="$(echo "$inst_line" | grep -oP '^[0-9]+')"
        edges="$(echo "$inst_line" | grep -oP '[0-9]+(?= edges)')"
    fi

    # Parse Objective line
    obj="" bound="" gap_pct="" time_s="" bb_nodes=""
    if obj_line="$(echo "$output" | grep -P '^Objective:')"; then
        obj="$(round3 "$(echo "$obj_line" | grep -oP 'Objective: \K[-0-9.e+]+')")"
        bound="$(round3 "$(echo "$obj_line" | grep -oP 'Bound: \K[-0-9.e+]+')")"
        gap_pct="$(round3 "$(echo "$obj_line" | grep -oP 'Gap: \K[0-9.e+-]+')")"
        time_s="$(round3 "$(echo "$obj_line" | grep -oP 'Time: \K[0-9.]+')")"
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
            stime="$(round3 "$(echo "$line" | grep -oP '[0-9.]+(?=s\s*$)')")" || stime=""
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
        warmstart_ub="$(round3 "$(echo "$ws_line" | awk '{print $3}')")"
        warmstart_time="$(round3 "$(echo "$ws_line" | awk '{gsub(/s$/,"",$5); print $5}')")"
    fi

    # LP iterations: total, plus the strong-branching / separation / heuristic
    # split, printed by cptp-solve under the "LP iterations" summary line.
    lp_iters="$(echo "$output" | grep -oP 'LP iterations\s+\K[0-9]+')" || lp_iters=""
    lp_iters_strongbr="$(echo "$output" | grep -oP '[0-9]+(?= \(strong br\.\))')" || lp_iters_strongbr=""
    lp_iters_sep="$(echo "$output" | grep -oP '[0-9]+(?= \(separation\))')" || lp_iters_sep=""
    lp_iters_heur="$(echo "$output" | grep -oP '[0-9]+(?= \(heuristics\))')" || lp_iters_heur=""

    # Parse per-phase times and fixing frequencies from the structured block.
    phases="$(echo "$output" | grep -P '^Phases:')" || phases=""
    prop="$(echo "$output" | grep -P '^Propagation:')" || prop=""
    rcfix="$(echo "$output" | grep -P '^RC-fixing:')" || rcfix=""
    sep_time="$(round3 "$(echo "$phases" | grep -oP 'separation \K[0-9.]+')")" || sep_time=""
    prop_time="$(round3 "$(echo "$phases" | grep -oP 'propagation \K[0-9.]+')")" || prop_time=""
    heur_time="$(round3 "$(echo "$phases" | grep -oP 'heuristic \K[0-9.]+')")" || heur_time=""
    rcfix_time="$(round3 "$(echo "$phases" | grep -oP 'rc-fixing \K[0-9.]+')")" || rcfix_time=""
    sweep_fixings="$(echo "$prop" | grep -oP 'sweep \K[0-9]+')" || sweep_fixings=""
    chain_fixings="$(echo "$prop" | grep -oP 'chain \K[0-9]+')" || chain_fixings=""
    node_fixings="$(echo "$prop" | grep -oP 'node \K[0-9]+')" || node_fixings=""
    rc_fix0="$(echo "$rcfix" | grep -oP 'fix0 \K[0-9]+')" || rc_fix0=""
    rc_fix1="$(echo "$rcfix" | grep -oP 'fix1 \K[0-9]+')" || rc_fix1=""

    # Write CSV row with timestamp and machine
    ROW_TIMESTAMP="$(date -Iseconds)"
    echo "${name},${nodes_graph},${edges},${obj},${bound},${gap_pct},${time_s},${TIME_LIMIT},${bb_nodes},${lp_iters},${lp_iters_strongbr},${lp_iters_sep},${lp_iters_heur},${total_cuts},${cut_rounds},${sec_stats},${rci_stats},${multistar_stats},${rglm_stats},${comb_stats},${warmstart_ub},${warmstart_time},${ROW_TIMESTAMP},\"${MACHINE}\",${sep_time},${prop_time},${heur_time},${rcfix_time},${sweep_fixings},${chain_fixings},${node_fixings},${rc_fix0},${rc_fix1}" >> "$OUTFILE"
    DONE["$name"]=1
    ran_count=$((ran_count + 1))

    # Progress line
    printf "%-35s %12s %10s %8s %7s\n" "$name" "${obj:-N/A}" "${time_s:+${time_s}s}" "${bb_nodes:-N/A}" "${total_cuts:-N/A}"
done

printf '%s\n' "$(printf '%.0s-' {1..75})"
# Summary: how much of the requested set is done, and whether to re-run.
total_targets=${#INSTANCE_STEMS[@]}
done_targets=0
for stem in "${INSTANCE_STEMS[@]}"; do
    [[ -n "${DONE[$stem]:-}" ]] && done_targets=$((done_targets + 1))
done
remaining=$((total_targets - done_targets))
echo "Ran ${ran_count} this invocation; skipped ${skipped_count} already-done."
echo "Completed ${done_targets}/${total_targets} requested instances; ${remaining} remaining."
if [[ "$stopped_early" -eq 1 || "$remaining" -gt 0 ]]; then
    echo "Re-run the same command to continue."
fi
echo "Done. Results written to $OUTFILE"
