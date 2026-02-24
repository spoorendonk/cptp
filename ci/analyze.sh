#!/usr/bin/env bash
# Main CI analysis entrypoint.
# Builds the project, captures warnings/errors, runs tests,
# then feeds failures to the local LLM for analysis.
# Outputs structured ISSUE_TITLE / ISSUE_BODY for the GitHub workflow.
set -euo pipefail

MODEL_PATH="${MODEL_PATH:-/models/model.gguf}"
LLAMA_PORT="${LLAMA_PORT:-8012}"
CONTEXT_SIZE="${CONTEXT_SIZE:-4096}"
MODE="${1:-auto}"  # auto | review

echo "=== rcspp LLM analysis (mode: $MODE) ==="

# ── Start llama.cpp server in background ─────────────────────────────
echo "Starting llama.cpp server..."
llama-server \
    --model "$MODEL_PATH" \
    --port "$LLAMA_PORT" \
    --ctx-size "$CONTEXT_SIZE" \
    --threads "$(nproc)" \
    --log-disable \
    &
LLAMA_PID=$!

for i in $(seq 1 30); do
    if curl -sf "http://localhost:${LLAMA_PORT}/health" > /dev/null 2>&1; then
        echo "llama.cpp server ready."
        break
    fi
    if [ "$i" -eq 30 ]; then
        echo "ERROR: llama.cpp server failed to start" >&2
        exit 1
    fi
    sleep 1
done

cleanup() {
    kill "$LLAMA_PID" 2>/dev/null || true
    wait "$LLAMA_PID" 2>/dev/null || true
}
trap cleanup EXIT

# ── Build ────────────────────────────────────────────────────────────
BUILD_LOG=$(mktemp)
echo "Building project..."
BUILD_OK=true
if ! cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON 2>&1 | tee "$BUILD_LOG"; then
    BUILD_OK=false
elif ! cmake --build build -j"$(nproc)" 2>&1 | tee -a "$BUILD_LOG"; then
    BUILD_OK=false
fi

WARNINGS=$(grep -iE '(warning|error):' "$BUILD_LOG" | head -100 || true)

# ── Test ─────────────────────────────────────────────────────────────
TEST_LOG=$(mktemp)
TEST_OK=true
if [ "$BUILD_OK" = true ]; then
    echo "Running algorithm tests..."
    if [ -f build/rcspp_algo_tests ]; then
        ./build/rcspp_algo_tests --reporter compact 2>&1 | tee "$TEST_LOG" || TEST_OK=false
    fi
    if [ -f build/rcspp_tests ]; then
        echo "Running integration tests..."
        ./build/rcspp_tests --reporter compact 2>&1 | tee -a "$TEST_LOG" || TEST_OK=false
    fi
fi

TEST_FAILURES=$(grep -E '(FAILED|ERROR)' "$TEST_LOG" | head -50 || true)

# ── Analyze with LLM ────────────────────────────────────────────────
ANALYSIS=$(mktemp)
ISSUE_TYPE=""

if [ "$MODE" = "review" ]; then
    ISSUE_TYPE="review"
    llm-review.sh review "" "$LLAMA_PORT" > "$ANALYSIS"

elif [ "$BUILD_OK" = false ]; then
    ISSUE_TYPE="build_failure"
    llm-review.sh build_failure "$WARNINGS" "$LLAMA_PORT" > "$ANALYSIS"

elif [ "$TEST_OK" = false ]; then
    ISSUE_TYPE="test_failure"
    llm-review.sh test_failure "$TEST_FAILURES" "$LLAMA_PORT" > "$ANALYSIS"

elif [ -n "$WARNINGS" ]; then
    ISSUE_TYPE="warnings"
    llm-review.sh warnings "$WARNINGS" "$LLAMA_PORT" > "$ANALYSIS"

else
    echo "Build and tests passed cleanly."
    if [ "$MODE" = "auto" ]; then
        echo "No issues found. Skipping LLM analysis."
        exit 0
    fi
    ISSUE_TYPE="review"
    llm-review.sh review "" "$LLAMA_PORT" > "$ANALYSIS"
fi

echo ""
echo "=== LLM Analysis ==="
cat "$ANALYSIS"

# ── Emit structured output for GitHub Actions ────────────────────────
COMMIT=$(cd /workspace && git rev-parse --short HEAD 2>/dev/null || echo "unknown")

case "$ISSUE_TYPE" in
    build_failure) TITLE="Build failure on ${COMMIT}" ;;
    test_failure)  TITLE="Test failure on ${COMMIT}" ;;
    warnings)      TITLE="Compiler warnings on ${COMMIT}" ;;
    review)        TITLE="Code review for ${COMMIT}" ;;
esac

echo "ISSUE_TITLE=${TITLE}"
echo "ISSUE_BODY<<EOF"
cat <<BODY
## ${TITLE}

**Commit:** \`${COMMIT}\`
**Mode:** ${MODE}
**Build:** $([ "$BUILD_OK" = true ] && echo "passed" || echo "FAILED")
**Tests:** $([ "$TEST_OK" = true ] && echo "passed" || echo "FAILED")

### Analysis

$(cat "$ANALYSIS")

---
_Automated analysis by Qwen2.5-Coder-1.5B (CPU) via llama.cpp_
BODY
echo "EOF"

if [ "$BUILD_OK" = false ] || [ "$TEST_OK" = false ]; then
    exit 1
fi
