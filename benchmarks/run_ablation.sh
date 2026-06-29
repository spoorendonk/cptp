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
#   ./benchmarks/run_ablation.sh [--time-limit N] [--max-hours H] [--resume] \
#                                [--configs "a_sec b_cap ..."] [instance_or_dir ...]
#
#   --time-limit N    Per-instance wall-clock limit in seconds (default: 3600).
#   --max-hours H     Overall wall-clock budget for THIS invocation in hours
#                     (float, default: 0 = unlimited). Run the sweep "H hours at a
#                     time": stops once H hours elapse and skips already-completed
#                     (config, instance) rows, so re-running the same command
#                     continues where it left off. Implies --resume.
#   --resume          Skip (config, instance) pairs already present in the CSV and
#                     run only the missing ones (no clock). Without --max-hours/
#                     --resume, existing rows for those pairs are replaced.
#   --configs "..."   Space-separated subset of config names (default: all below).
#   instance_or_dir   Instance file(s)/dir(s). Default: benchmarks/instances/{spprclib,roberti}.
#
# Output: benchmarks/ablation.csv, plus benchmarks/logs/ablation/<config>/<instance>.log.gz
#   (full per-run stdout, gitignored). Default mode replaces rows for the requested
#   (config, instance) pairs; --resume/--max-hours skip already-completed pairs. If
#   the CSV header is out of date, the old file is archived to ablation.csv.bak first.
#
# Config table (name : extra flags passed to cptp-solve):
#   a_sec        SEC only                         --enable_rci false --enable_multistar false
#   b_cap        + capacity cuts (default set)    (no extra flags: SEC+RCI+Multistar)
#   b_cap_nofix  b_cap, fixing/propagation OFF    --rc_fixing off --edge_elimination false
#                                                 --bounds_propagation false (no comb/rglm)
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
MAX_HOURS=0
RESUME=0
CONFIGS_REQUESTED=""
TARGETS=()

config_flags() {
    case "$1" in
        a_sec)      echo "--enable_rci false --enable_multistar false" ;;
        b_cap)      echo "" ;;
        b_cap_nofix) echo "--rc_fixing off --edge_elimination false --bounds_propagation false" ;;
        c_combrglm) echo "--enable_comb true --enable_rglm true" ;;
        d_nofix)    echo "--enable_comb true --enable_rglm true --rc_fixing off --edge_elimination false --bounds_propagation false" ;;
        d_fix)      echo "--enable_comb true --enable_rglm true --rc_fixing adaptive" ;;
        spi_off)    echo "--all_pairs_bounds true" ;;
        spi_on)     echo "--all_pairs_bounds true --enable_spi true" ;;
        *)          echo "__UNKNOWN__" ;;
    esac
}
ALL_CONFIGS="a_sec b_cap b_cap_nofix c_combrglm d_nofix d_fix spi_off spi_on"

usage() { sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-limit) TIME_LIMIT="$2"; shift 2 ;;
        --max-hours)  MAX_HOURS="$2"; shift 2 ;;
        --resume)     RESUME=1; shift ;;
        --configs)    CONFIGS_REQUESTED="$2"; shift 2 ;;
        --help|-h)    usage 0 ;;
        -*)           echo "Unknown option: $1" >&2; usage 1 ;;
        *)            TARGETS+=("$1"); shift ;;
    esac
done

# Wall-clock budget (seconds) for this invocation; 0 = unlimited. A budget
# implies resume/skip-completed so repeated "H hours at a time" runs progress.
# Validate up front so a non-numeric/negative value errors instead of silently
# parsing to 0 (= unlimited) and running unbounded.
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
LOGROOT="$SCRIPT_DIR/logs/ablation"
# New phase/fixing columns are appended after `machine` (append-only). sep_time
# is the separation-phase wall-clock (one span per round around the parallel
# separators); prop_time includes rcfix_time (RC fixing runs inside the
# propagator); heur_time is the LP-guided B&C callback only.
HEADER="config,instance,nodes_graph,edges,obj,bound,gap_pct,time_s,time_limit,bb_nodes,total_cuts,cut_rounds,sec_cuts,sec_rounds,sec_time,rci_cuts,rci_rounds,rci_time,multistar_cuts,multistar_rounds,multistar_time,rglm_cuts,rglm_rounds,rglm_time,comb_cuts,comb_rounds,comb_time,spi_cuts,spi_rounds,spi_time,timestamp,machine,sep_time,prop_time,heur_time,rcfix_time,sweep_fixings,chain_fixings,node_fixings,rc_fix0,rc_fix1"
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

MACHINE="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g') $(nproc)T"

# Drop any existing CSV row for a (config,instance) pair (exact match on fields
# 1-2), so a fresh run replaces it rather than appending a duplicate. Called
# per-pair right before a run, so a budget-stop leaves not-yet-run rows (incl.
# widened/padded ones) intact. Atomic via temp + mv.
drop_row() {
    local cfg="$1" inst="$2" tmp
    tmp="$(mktemp)"
    awk -F, -v c="$cfg" -v i="$inst" 'NR==1 || !($1==c && $2==i)' "$OUTFILE" > "$tmp"
    mv "$tmp" "$OUTFILE"
}

# Completed (config,instance) set for resume/skip.
# A row counts as completed only if its new columns are populated (last field,
# rc_fix1, non-empty); rows widened from an old schema have empty new fields, so
# resume re-runs them to fill the new columns instead of skipping them.
declare -A DONE
while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    [[ -n "${line##*,}" ]] || continue
    DONE["$(echo "$line" | cut -d, -f1,2)"]=1
done < <(tail -n +2 "$OUTFILE")

RUN_START="$(date +%s)"
ran_count=0
skipped_count=0
stopped_early=0

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

        # Stop before starting a new run once the wall budget is spent (an
        # already-running instance still finishes; overshoot <= one --time-limit).
        if [[ "$BUDGET_SEC" -gt 0 && $(( $(date +%s) - RUN_START )) -ge "$BUDGET_SEC" ]]; then
            echo "Wall budget ${MAX_HOURS}h reached — stopping; re-run to continue."
            stopped_early=1
            break 2
        fi

        # Resume/skip: leave already-completed (config,instance) rows untouched.
        if [[ "$RESUME_SKIP" -eq 1 && -n "${DONE[$cfg,$name]:-}" ]]; then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        # Replace any existing row for this pair (stale, or widened/padded).
        drop_row "$cfg" "$name"

        # Run solver; persist full stdout, parse from the saved file (so CSV and
        # log agree), then gzip. parse_sep() reads the $OUTPUT global.
        logdir="$LOGROOT/$cfg"
        mkdir -p "$logdir"
        logfile="$logdir/$name.log"
        # shellcheck disable=SC2086
        "$SOLVER" "$inst" --time_limit "$TIME_LIMIT" $flags > "$logfile" 2>&1 || true
        OUTPUT="$(cat "$logfile")"
        gzip -f "$logfile"

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

        # Per-phase times and fixing frequencies from the structured block.
        phases="$(echo "$OUTPUT" | grep -P '^Phases:')" || phases=""
        prop="$(echo "$OUTPUT" | grep -P '^Propagation:')" || prop=""
        rcfix="$(echo "$OUTPUT" | grep -P '^RC-fixing:')" || rcfix=""
        sep_time="$(round3 "$(echo "$phases" | grep -oP 'separation \K[0-9.]+')")" || sep_time=""
        prop_time="$(round3 "$(echo "$phases" | grep -oP 'propagation \K[0-9.]+')")" || prop_time=""
        heur_time="$(round3 "$(echo "$phases" | grep -oP 'heuristic \K[0-9.]+')")" || heur_time=""
        rcfix_time="$(round3 "$(echo "$phases" | grep -oP 'rc-fixing \K[0-9.]+')")" || rcfix_time=""
        sweep_fixings="$(echo "$prop" | grep -oP 'sweep \K[0-9]+')" || sweep_fixings=""
        chain_fixings="$(echo "$prop" | grep -oP 'chain \K[0-9]+')" || chain_fixings=""
        node_fixings="$(echo "$prop" | grep -oP 'node \K[0-9]+')" || node_fixings=""
        rc_fix0="$(echo "$rcfix" | grep -oP 'fix0 \K[0-9]+')" || rc_fix0=""
        rc_fix1="$(echo "$rcfix" | grep -oP 'fix1 \K[0-9]+')" || rc_fix1=""
        ts="$(date -Iseconds)"

        echo "${cfg},${name},${nodes_graph},${edges},${obj},${bound},${gap_pct},${time_s},${TIME_LIMIT},${bb_nodes},${total_cuts},${cut_rounds},${sec_stats},${rci_stats},${multistar_stats},${rglm_stats},${comb_stats},${spi_stats},${ts},\"${MACHINE}\",${sep_time},${prop_time},${heur_time},${rcfix_time},${sweep_fixings},${chain_fixings},${node_fixings},${rc_fix0},${rc_fix1}" >> "$OUTFILE"
        DONE["$cfg,$name"]=1
        ran_count=$((ran_count + 1))
        printf "%-10s %-30s %12s %9s %8s\n" "$cfg" "$name" "${obj:-NA}" "${time_s:-NA}" "${bb_nodes:-NA}"
    done
done

# Summary across the requested (config,instance) grid.
total_targets=0; done_targets=0
for cfg in $CONFIGS; do
    for inst in "${INSTANCES[@]}"; do
        total_targets=$((total_targets + 1))
        [[ -n "${DONE[$cfg,$(basename "${inst%.*}")]:-}" ]] && done_targets=$((done_targets + 1))
    done
done
remaining=$((total_targets - done_targets))
echo "Ran ${ran_count} this invocation; skipped ${skipped_count} already-done."
echo "Completed ${done_targets}/${total_targets} requested (config,instance) pairs; ${remaining} remaining."
if [[ "$stopped_early" -eq 1 || "$remaining" -gt 0 ]]; then
    echo "Re-run the same command to continue."
fi
echo "Done. Rows written to $OUTFILE"
