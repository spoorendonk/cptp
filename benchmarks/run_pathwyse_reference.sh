#!/usr/bin/env bash
# run_pathwyse_reference.sh — Generate PathWyse reference results for benchmarking.
#
# Runs PathWyse on all instances (spprclib + roberti) via compare_pathwyse.py --csv.
# Instance conversion (.sppcc/.vrp → PathWyse native format) is handled internally
# by compare_pathwyse.py.
#
# Usage:
#   ./benchmarks/run_pathwyse_reference.sh [--time-limit N] [--algo A] [--out F]
#                                          [instance_or_dir ...]
#
# Arguments:
#   --time-limit N       Per-instance time limit in seconds (default: 3600)
#   --algo A             PWDefault (default) or PWBucket (PathWyse's suggested
#                        settings for cyclic networks)
#   --out F              Output CSV (default: benchmarks/pathwyse-<algo>.csv)
#   instance_or_dir      Instance files or directories of instances, e.g.
#                        benchmarks/instances/multicap/spprclib/m2 for the
#                        multi-capacity sets. Default: both
#                        benchmarks/instances/spprclib and
#                        benchmarks/instances/roberti.
#
# Output:
#   The --out CSV, with algo, timestamp, machine and PathWyse commit columns.
#   Existing rows for re-run instances are replaced; other rows are preserved.
#   Each target's rows replace the old ones only once that target has run to
#   the end, so a target that fails or is interrupted keeps its old rows.
#   benchmarks/pathwyse.csv holds the paper's rows (an unpinned PathWyse, see
#   setup_pathwyse.sh) in the older column layout and is not written to.
#
# Prerequisites:
#   PathWyse must be built: run ./benchmarks/setup_pathwyse.sh first. Only the
#   pinned build runs: a binary whose provenance stamp does not prove it was
#   built at setup_pathwyse.sh's PATHWYSE_COMMIT is refused (PATHWYSE_BIN may
#   point elsewhere, at a binary with such a stamp beside it).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TIME_LIMIT=3600
ALGO=PWDefault
OUTFILE=""
TARGETS=()

usage() {
    sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-limit)
            TIME_LIMIT="$2"; shift 2 ;;
        --algo)
            ALGO="$2"; shift 2 ;;
        --out)
            OUTFILE="$2"; shift 2 ;;
        --help|-h)
            usage 0 ;;
        -*)
            echo "Unknown option: $1" >&2; usage 1 ;;
        *)
            TARGETS+=("$1"); shift ;;
    esac
done

# Default targets: both instance directories
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=(
        "$SCRIPT_DIR/instances/spprclib"
        "$SCRIPT_DIR/instances/roberti"
    )
fi

if [[ -z "$OUTFILE" ]]; then
    OUTFILE="$SCRIPT_DIR/pathwyse-$ALGO.csv"
fi

# Refuse a file in another column layout (e.g. pathwyse.csv) up front, as
# compare_pathwyse.py would for each target; it writes this header.
EXPECTED_HEADER="instance,set,algo,obj,time_s,status,timestamp,machine,pathwyse,note"
if [[ -f "$OUTFILE" ]]; then
    have_header="$(grep -v '^#' "$OUTFILE" | grep -m1 . || true)"
    have_header="${have_header%$'\r'}"  # Python's csv writer ends lines with \r\n
    if [[ -n "$have_header" && "$have_header" != "$EXPECTED_HEADER" ]]; then
        echo "Error: $OUTFILE has header '$have_header', expected '$EXPECTED_HEADER'" >&2
        exit 1
    fi
fi

# Detect machine info
if command -v lscpu >/dev/null; then
    CPU="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g')"
else
    CPU="$(sysctl -n machdep.cpu.brand_string)"
fi
if command -v nproc >/dev/null; then
    NCPU="$(nproc)"
else
    NCPU="$(sysctl -n hw.logicalcpu)"
fi
MACHINE="$CPU ${NCPU}T"

# The pin setup_pathwyse.sh builds; compare_pathwyse.py refuses any other build.
PIN="$(sed -n 's/^PATHWYSE_COMMIT="\([0-9a-f]\{40\}\)"$/\1/p' "$SCRIPT_DIR/setup_pathwyse.sh")"
if [[ -z "$PIN" ]]; then
    echo "Error: no PATHWYSE_COMMIT pin found in $SCRIPT_DIR/setup_pathwyse.sh" >&2
    exit 1
fi

echo "PathWyse reference run — algorithm: $ALGO, time limit: ${TIME_LIMIT}s"
echo "PathWyse pin: $PIN"
echo "Output: $OUTFILE"
echo

# compare_pathwyse.py replaces a target's old rows only after the target has
# run, so nothing is removed up front. A target with a mismatch exits 1 after
# writing its rows: run the rest anyway and fail at the end. A refusal before
# any run (exit 2: an unproven binary, a wrong header, a missing target)
# stops the whole run, with the rows of the targets before it written.
STATUS=0
for target in "${TARGETS[@]}"; do
    echo "--- Running: $target ---"
    rc=0
    python3 "$SCRIPT_DIR/compare_pathwyse.py" "$target" \
        --time-limit "$TIME_LIMIT" \
        --algo "$ALGO" \
        --csv "$OUTFILE" \
        --machine "$MACHINE" \
        --require-pathwyse "$PIN" || rc=$?
    if [[ $rc -eq 2 ]]; then
        echo "Error: compare_pathwyse.py refused $target; $OUTFILE is unchanged by it" >&2
        exit 2
    elif [[ $rc -ne 0 ]]; then
        STATUS=1
    fi
    echo
done

echo "Done. Results written to $OUTFILE"
exit "$STATUS"
