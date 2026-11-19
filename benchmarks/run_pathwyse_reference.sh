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
#   benchmarks/pathwyse.csv holds the paper's rows (an unpinned PathWyse, see
#   setup_pathwyse.sh) in the older column layout and is not written to.
#
# Prerequisites:
#   PathWyse must be built: run ./benchmarks/setup_pathwyse.sh first.

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

# Refuse a file in another column layout (e.g. pathwyse.csv) before the
# row replacement below rewrites it; compare_pathwyse.py writes this header.
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

# Collect instance stems to be (re-)run for dedup
collect_stems() {
    local target="$1"
    if [[ -f "$target" ]]; then
        basename "${target%.*}"
    elif [[ -d "$target" ]]; then
        find "$target" -maxdepth 1 \( -name '*.sppcc' -o -name '*.vrp' \) -exec basename {} \; | sed 's/\.[^.]*$//' | sort
    fi
}

INSTANCE_STEMS=()
for target in "${TARGETS[@]}"; do
    while IFS= read -r stem; do
        INSTANCE_STEMS+=("$stem")
    done < <(collect_stems "$target")
done

# Remove existing rows for instances we're about to run
if [[ -f "$OUTFILE" ]]; then
    TMPFILE="$(mktemp)"
    trap 'rm -f "$TMPFILE"' EXIT
    head -1 "$OUTFILE" > "$TMPFILE"
    tail -n +2 "$OUTFILE" | while IFS= read -r line; do
        # Replace only this algorithm's rows: a file may hold several
        IFS=, read -r row_instance _ row_algo _ <<< "$line"
        skip=false
        for stem in "${INSTANCE_STEMS[@]}"; do
            if [[ "$row_instance" == "$stem" && "$row_algo" == "$ALGO" ]]; then
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
fi

echo "PathWyse reference run — algorithm: $ALGO, time limit: ${TIME_LIMIT}s"
echo "Output: $OUTFILE"
echo

# A target with a mismatch exits 1 after writing its rows; run the rest
# anyway (their old rows are already gone) and fail at the end.
STATUS=0
for target in "${TARGETS[@]}"; do
    echo "--- Running: $target ---"
    python3 "$SCRIPT_DIR/compare_pathwyse.py" "$target" \
        --time-limit "$TIME_LIMIT" \
        --algo "$ALGO" \
        --csv "$OUTFILE" \
        --machine "$MACHINE" || STATUS=1
    echo
done

echo "Done. Results written to $OUTFILE"
exit "$STATUS"
