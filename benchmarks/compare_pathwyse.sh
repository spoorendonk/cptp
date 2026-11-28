#!/usr/bin/env bash
# compare_pathwyse.sh — Compare PathWyse and cptp benchmark results.
#
# Reads benchmarks/pathwyse.csv and benchmarks/cptp.csv and prints a comparison
# table with objective match, time speedup, and summary statistics.
#
# Usage:
#   ./benchmarks/compare_pathwyse.sh [--pathwyse <file>] [--cptp <file>]
#
# Arguments:
#   --pathwyse <file>    PathWyse results CSV (default: benchmarks/pathwyse.csv,
#                        the paper's rows; run_pathwyse_reference.sh writes
#                        benchmarks/pathwyse-<algo>.csv)
#   --cptp <file>        cptp results CSV (default: benchmarks/cptp.csv)
#
# Output:
#   Comparison table with columns:
#     instance, pw_obj, pw_time, cptp_obj, cptp_time, speedup, obj_match
#   A pair is compared only when both rows prove optimality (PathWyse status
#   Optimal; cptp gap_pct <= 0.01, HiGHS's default relative gap as a percent,
#   or obj equal to bound before time_limit): OK or FAIL. Otherwise obj_match says why not:
#   TLIM (time limit), ERR (failed run), UNPR (ended unproven), ? (no row).
#   Summary: match/mismatch counts, proven-optimal counts per solver, and the
#   geometric mean speedup over the OK pairs only. Exits 1 on a FAIL.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PW_FILE="" CPTP_FILE=""

usage() {
    sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pathwyse) PW_FILE="$2"; shift 2 ;;
        --cptp)     CPTP_FILE="$2"; shift 2 ;;
        --help|-h)  usage 0 ;;
        *)          echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

# Default to committed CSV files
if [[ -z "$PW_FILE" ]]; then
    PW_FILE="$SCRIPT_DIR/pathwyse.csv"
    if [[ ! -f "$PW_FILE" ]]; then
        echo "Error: $PW_FILE not found" >&2
        echo "It holds the paper's rows; pass --pathwyse benchmarks/pathwyse-<algo>.csv" >&2
        echo "for a run_pathwyse_reference.sh run." >&2
        exit 1
    fi
fi
if [[ -z "$CPTP_FILE" ]]; then
    CPTP_FILE="$SCRIPT_DIR/cptp.csv"
    if [[ ! -f "$CPTP_FILE" ]]; then
        echo "Error: $CPTP_FILE not found" >&2
        echo "Run ./benchmarks/run_benchmarks.sh first." >&2
        exit 1
    fi
fi

echo "PathWyse: $PW_FILE"
echo "cptp:     $CPTP_FILE"
echo

python3 - "$PW_FILE" "$CPTP_FILE" << 'PYEOF'
import csv
import math
import sys

pw_file, cptp_file = sys.argv[1], sys.argv[2]

def read_csv(path):
    """Read CSV, skipping comment lines."""
    rows = {}
    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]
    reader = csv.DictReader(lines)
    algos = set()
    for row in reader:
        rows[row['instance']] = row
        if row.get('algo'):
            algos.add(row['algo'])
    # Rows are keyed by instance: a file holding several PathWyse algorithms
    # (run_pathwyse_reference.sh --out shared by two --algo runs) would mix them.
    if len(algos) > 1:
        sys.exit(f"Error: {path} holds rows of several algorithms {sorted(algos)}; "
                 "compare one algorithm's file at a time")
    return rows

pw = read_csv(pw_file)
cptp = read_csv(cptp_file)
if any('gap_pct' not in row for row in cptp.values()):
    sys.exit(f"Error: {cptp_file} has no gap_pct column, so no cptp row proves optimality")

# Merge on instance name
all_instances = sorted(set(pw.keys()) | set(cptp.keys()))

header = f"{'Instance':<35} {'PW Obj':>12} {'PW Time':>10} {'CPTP Obj':>12} {'CPTP Time':>10} {'Speedup':>8} {'Match':>6}"
print(header)
print('-' * len(header))

def number(row, key):
    value = row.get(key, '')
    return float(value) if value else None

def pw_verdict(row):
    """'' if the PathWyse row proves optimality, else why not.

    Only status Optimal proves it (a failed run's row keeps no objective,
    and one that ended Feasible or Super optimal proves nothing)."""
    status = row.get('status', '')
    if status == 'TIMEOUT':
        return 'TLIM'
    if not status or status.startswith('ERROR') or number(row, 'obj') is None:
        return 'ERR'
    return '' if status == 'Optimal' else 'UNPR'

def cptp_verdict(row):
    """'' if the cptp row proves optimality, else why not.

    Proven is an objective at a gap_pct (a percent, #187) within HiGHS's
    default relative gap, 1e-4, or one equal to its bound as recorded in a
    run that ended before its time limit: HiGHS also stops at its absolute
    gap (1e-6), which on an objective near 0 is a relative gap above 0.01%.
    (obj and bound are rounded to 3 decimals, so equal alone proves
    nothing.)"""
    gap = number(row, 'gap_pct')
    obj = number(row, 'obj')
    if obj is None or gap is None:
        return 'ERR'
    time_s, time_limit = number(row, 'time_s'), number(row, 'time_limit')
    in_time = time_s is not None and time_limit is not None and time_s < time_limit
    at_bound = in_time and obj == number(row, 'bound')
    return '' if gap <= 0.01 or at_bound else 'TLIM'

matches = mismatches = unproven = missing = 0
pw_solved = cptp_solved = 0
log_speedups = []

for name in all_instances:
    p = pw.get(name)
    c = cptp.get(name)
    p_why = pw_verdict(p) if p is not None else '?'
    c_why = cptp_verdict(c) if c is not None else '?'
    pw_solved += p_why == ''
    cptp_solved += c_why == ''
    p, c = p or {}, c or {}

    # Only a pair that both prove optimal is compared, and only an agreeing
    # pair counts toward the speedup: a failed, timed-out or unproven row is
    # neither a match nor a speedup.
    speedup_str = ''
    if p_why == '?' or c_why == '?':
        match_str = '?'
        missing += 1
    elif p_why or c_why:
        match_str = p_why or c_why
        unproven += 1
    elif abs(number(p, 'obj') - number(c, 'obj')) < 0.01:
        match_str = 'OK'
        matches += 1
        pw_time, cptp_time = number(p, 'time_s'), number(c, 'time_s')
        if pw_time and cptp_time:
            speedup = pw_time / cptp_time
            speedup_str = f'{speedup:.2f}x'
            log_speedups.append(math.log(speedup))
    else:
        match_str = 'FAIL'
        mismatches += 1

    print(f"{name:<35} {p.get('obj', ''):>12} {p.get('time_s', ''):>10} {c.get('obj', ''):>12} {c.get('time_s', ''):>10} {speedup_str:>8} {match_str:>6}")

print('-' * len(header))

# Summary
total = matches + mismatches + unproven + missing
print(f"\nInstances: {total}  Matches: {matches}  Mismatches: {mismatches}  "
      f"Not both proven: {unproven}  Missing: {missing}")
print(f"Proven optimal: PathWyse {pw_solved}, cptp {cptp_solved}")
if log_speedups:
    geo_mean = math.exp(sum(log_speedups) / len(log_speedups))
    print(f"Geometric mean speedup (PW/cptp): {geo_mean:.2f}x  (over {len(log_speedups)} matching instances)")
if mismatches:
    sys.exit(1)
PYEOF
