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
#   --pathwyse <file>    PathWyse results CSV (default: benchmarks/pathwyse.csv)
#   --cptp <file>        cptp results CSV (default: benchmarks/cptp.csv)
#
# Output:
#   Comparison table with columns:
#     instance, pw_obj, pw_time, cptp_obj, cptp_time, speedup, obj_match
#   Summary: geometric mean speedup, match/mismatch counts.

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
        echo "Run ./benchmarks/run_pathwyse_reference.sh first." >&2
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
    for row in reader:
        rows[row['instance']] = row
    return rows

pw = read_csv(pw_file)
cptp = read_csv(cptp_file)

# Merge on instance name
all_instances = sorted(set(pw.keys()) | set(cptp.keys()))

header = f"{'Instance':<35} {'PW Obj':>12} {'PW Time':>10} {'CPTP Obj':>12} {'CPTP Time':>10} {'Speedup':>8} {'Match':>6}"
print(header)
print('-' * len(header))

matches = mismatches = missing = 0
log_speedups = []

for name in all_instances:
    p = pw.get(name, {})
    c = cptp.get(name, {})

    pw_obj_str = p.get('obj', '')
    pw_time_str = p.get('time_s', '')
    cptp_obj_str = c.get('obj', '')
    cptp_time_str = c.get('time_s', '')

    pw_obj = float(pw_obj_str) if pw_obj_str else None
    pw_time = float(pw_time_str) if pw_time_str else None
    cptp_obj = float(cptp_obj_str) if cptp_obj_str else None
    cptp_time = float(cptp_time_str) if cptp_time_str else None

    # Objective match
    match_str = '?'
    if pw_obj is not None and cptp_obj is not None:
        if abs(pw_obj - cptp_obj) < 0.01:
            match_str = 'OK'
            matches += 1
        else:
            match_str = 'FAIL'
            mismatches += 1
    else:
        missing += 1

    # Speedup (PathWyse time / cptp time)
    speedup_str = ''
    if pw_time and cptp_time and cptp_time > 0:
        speedup = pw_time / cptp_time
        speedup_str = f'{speedup:.2f}x'
        if pw_time > 0:
            log_speedups.append(math.log(speedup))

    print(f"{name:<35} {pw_obj_str:>12} {pw_time_str:>10} {cptp_obj_str:>12} {cptp_time_str:>10} {speedup_str:>8} {match_str:>6}")

print('-' * len(header))

# Summary
total = matches + mismatches + missing
print(f"\nInstances: {total}  Matches: {matches}  Mismatches: {mismatches}  Missing: {missing}")
if log_speedups:
    geo_mean = math.exp(sum(log_speedups) / len(log_speedups))
    print(f"Geometric mean speedup (PW/cptp): {geo_mean:.2f}x  (over {len(log_speedups)} comparable instances)")
PYEOF
