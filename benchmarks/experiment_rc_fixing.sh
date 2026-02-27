#!/usr/bin/env bash
set -euo pipefail

# RC fixing policy comparison for Step 0.5.
# Output: benchmarks/rc_fixing_study.csv

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT_DIR/build/rcspp-solve"
OUT_CSV="$ROOT_DIR/benchmarks/rc_fixing_study.csv"

if [[ ! -x "$BIN" ]]; then
  echo "error: binary not found: $BIN" >&2
  echo "build first: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j" >&2
  exit 1
fi

INSTANCES=(
  "benchmarks/instances/spprclib/B-n45-k6-54.sppcc"
  "benchmarks/instances/spprclib/P-n50-k7-92.sppcc"
  "benchmarks/instances/spprclib/A-n63-k9-157.sppcc"
)

STRATEGIES=(off root_only on_ub_improvement adaptive)
DSSR=(false true)
TIME_LIMIT="${1:-60}"

printf "instance,dssr_async,rc_fixing,time_sec,nodes,obj,bound,gap_pct,rc_fix0,rc_fix1,rc_label_runs,rc_callback_runs,rc_time_sec\n" > "$OUT_CSV"

for inst in "${INSTANCES[@]}"; do
  for dssr in "${DSSR[@]}"; do
    for strat in "${STRATEGIES[@]}"; do
      echo "[run] $inst dssr_async=$dssr rc_fixing=$strat"
      out="$($BIN "$ROOT_DIR/$inst" \
        --time_limit "$TIME_LIMIT" \
        --output_flag false \
        --rc_fixing "$strat" \
        --dssr_async "$dssr" 2>&1)"

      summary_line="$(echo "$out" | rg "^Objective:" -n | tail -n1 | sed 's/^[0-9]*://' || true)"
      rc_line="$(echo "$out" | rg "^RC fixing:" -n | tail -n1 | sed 's/^[0-9]*://' || true)"

      obj="$(echo "$summary_line" | awk '{print $2}')"
      bound="$(echo "$summary_line" | awk '{print $4}')"
      gap="$(echo "$summary_line" | awk '{print $6}' | tr -d '%')"
      time_s="$(echo "$summary_line" | awk '{print $8}' | tr -d 's')"
      nodes="$(echo "$summary_line" | awk '{print $10}')"

      rc_fix0="0"; rc_fix1="0"; rc_runs="0"; rc_cb="0"; rc_t="0"
      if [[ -n "$rc_line" ]]; then
        rc_fix0="$(echo "$rc_line" | awk '{print $3}')"
        rc_fix1="$(echo "$rc_line" | awk '{print $8}')"
        rc_runs="$(echo "$rc_line" | awk '{print $13}' | tr -d ',')"
        rc_cb="$(echo "$rc_line" | awk '{print $16}' | tr -d ',')"
        rc_t="$(echo "$rc_line" | awk '{print $19}' | tr -d 's')"
      fi

      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$(basename "$inst")" "$dssr" "$strat" "$time_s" "$nodes" "$obj" "$bound" "$gap" \
        "$rc_fix0" "$rc_fix1" "$rc_runs" "$rc_cb" "$rc_t" >> "$OUT_CSV"
    done
  done
 done

echo "Wrote $OUT_CSV"
