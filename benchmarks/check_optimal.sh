#!/usr/bin/env bash
# check_optimal.sh — Check cptp benchmark results against known optimal values.
#
# Reads results from benchmarks/cptp.csv, looks up reference optima from
# benchmarks/instances/{spprclib,roberti}/optimal.csv, and reports per instance.
#
# Usage:
#   ./benchmarks/check_optimal.sh [--csv FILE] [PATH...]
#
# Arguments:
#   PATH           Instance file or directory — filter CSV to these instances.
#                  If omitted, all CSV rows are checked.
#   --csv FILE     Results CSV (default: benchmarks/cptp.csv).
#
# Verdict:
#   PASS   abs(obj - optimal) < 0.001
#   FAIL   wrong obj, or hit time limit (gap_pct > 0)
#   SKIP   no reference optimal available
#
# Output:
#   Per-instance table and PASS/FAIL/SKIP summary counts.
set -euo pipefail

usage() {
  sed -n '2,/^$/s/^# \?//p' "${BASH_SOURCE[0]}"
  exit "${1:-0}"
}

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
CSV="${SCRIPTDIR}/cptp.csv"
PATHS=()
EPS=0.001

while [[ $# -gt 0 ]]; do
  case "$1" in
    --csv)  CSV="$2"; shift 2 ;;
    -h|--help) usage ;;
    *)      PATHS+=("$1"); shift ;;
  esac
done

if [[ ! -f "$CSV" ]]; then
  echo "CSV not found: $CSV" >&2
  exit 1
fi

# ── Collect instance stems to filter ──
FILTER_STEMS=()
if [[ ${#PATHS[@]} -gt 0 ]]; then
  for p in "${PATHS[@]}"; do
    if [[ -f "$p" ]]; then
      stem="$(basename "$p")"
      FILTER_STEMS+=("${stem%.*}")
    elif [[ -d "$p" ]]; then
      while IFS= read -r -d '' f; do
        stem="$(basename "$f")"
        FILTER_STEMS+=("${stem%.*}")
      done < <(find "$p" -maxdepth 1 -type f \( -name '*.sppcc' -o -name '*.vrp' \) -print0)
    fi
  done
fi

# ── Load optimal values ──
declare -A OPTIMAL
for setdir in "$SCRIPTDIR"/instances/spprclib "$SCRIPTDIR"/instances/roberti; do
  optfile="$setdir/optimal.csv"
  [[ -f "$optfile" ]] || continue
  while IFS=, read -r inst opt; do
    [[ "$inst" == "instance" || "$inst" == \#* ]] && continue
    OPTIMAL["$inst"]="$opt"
  done < "$optfile"
done

# ── Check results ──
PASS=0 FAIL=0 SKIP=0 TOTAL=0

printf "%-30s  %10s  %10s  %s\n" "Instance" "Obj" "Optimal" "Result"
printf '%.0s-' {1..65}; echo

# Read CSV, skip header
while IFS=, read -r inst _ _ obj _ gap_pct _; do
  # Apply instance filter
  if [[ ${#FILTER_STEMS[@]} -gt 0 ]]; then
    match=0
    for s in "${FILTER_STEMS[@]}"; do
      if [[ "$inst" == "$s" ]]; then match=1; break; fi
    done
    [[ $match -eq 0 ]] && continue
  fi

  TOTAL=$((TOTAL + 1))

  # Time-limit check: gap_pct > 0 means not proven optimal
  is_tlim="$(awk -v g="$gap_pct" 'BEGIN{print (g+0 > 0) ? 1 : 0}')"
  if [[ "$is_tlim" -eq 1 ]]; then
    opt="${OPTIMAL[$inst]:-}"
    printf "%-30s  %10s  %10s  FAIL (TLIM)\n" "$inst" "$obj" "${opt:--}"
    FAIL=$((FAIL + 1))
    continue
  fi

  # Look up optimal
  opt="${OPTIMAL[$inst]:-}"
  if [[ -z "$opt" ]]; then
    printf "%-30s  %10s  %10s  SKIP\n" "$inst" "$obj" "-"
    SKIP=$((SKIP + 1))
    continue
  fi

  # Compare: abs(obj - optimal) < eps
  result="$(awk -v c="$obj" -v o="$opt" -v e="$EPS" 'BEGIN{ d=c-o; if(d<0)d=-d; print (d<e) ? "PASS" : "FAIL" }')"
  printf "%-30s  %10s  %10s  %s\n" "$inst" "$obj" "$opt" "$result"
  if [[ "$result" == "PASS" ]]; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
  fi
done < <(tail -n +2 "$CSV")

echo
echo "Summary: PASS=$PASS  FAIL=$FAIL  SKIP=$SKIP  TOTAL=$TOTAL"

[[ $FAIL -eq 0 ]]
