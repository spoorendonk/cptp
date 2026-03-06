#!/usr/bin/env bash
# run_pathwyse_reference.sh — Generate PathWyse reference results for benchmarking.
#
# Runs PathWyse on all instances (spprclib + roberti) via compare_pathwyse.py --csv.
# Instance conversion (.sppcc/.vrp → PathWyse native format) is handled internally
# by compare_pathwyse.py.
#
# Usage:
#   ./benchmarks/run_pathwyse_reference.sh [--time-limit N] [instance_or_dir]
#
# Arguments:
#   --time-limit N       Per-instance time limit in seconds (default: 3600)
#   instance_or_dir      Single instance file or directory of instances.
#                        Default: both benchmarks/instances/spprclib and
#                        benchmarks/instances/roberti.
#
# Output:
#   benchmarks/pathwyse.csv — committed CSV with timestamp and machine columns.
#   Existing rows for re-run instances are replaced; other rows are preserved.
#
# Prerequisites:
#   PathWyse must be built: run ./benchmarks/setup_pathwyse.sh first.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# Default targets: both instance directories
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=(
        "$SCRIPT_DIR/instances/spprclib"
        "$SCRIPT_DIR/instances/roberti"
    )
fi

OUTFILE="$SCRIPT_DIR/pathwyse.csv"

# Detect machine info
MACHINE="$(lscpu | grep 'Model name' | sed 's/.*: *//' | sed 's/  */ /g') $(nproc)T"

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
fi

echo "PathWyse reference run — time limit: ${TIME_LIMIT}s"
echo "Output: $OUTFILE"
echo

for target in "${TARGETS[@]}"; do
    echo "--- Running: $target ---"
    python3 "$SCRIPT_DIR/compare_pathwyse.py" "$target" \
        --time-limit "$TIME_LIMIT" \
        --csv "$OUTFILE" \
        --machine "$MACHINE"
    echo
done

echo "Done. Results written to $OUTFILE"
