#!/usr/bin/env bash
# Sends prompts to the local llama.cpp server for code analysis.
# Usage: llm-review.sh <mode> <context> <port>
set -euo pipefail

MODE="${1:-review}"
CONTEXT="${2:-}"
PORT="${3:-8012}"
API="http://localhost:${PORT}/v1/chat/completions"

query_llm() {
    local system_msg="$1"
    local user_msg="$2"

    # Escape for JSON
    system_msg=$(printf '%s' "$system_msg" | jq -Rsa .)
    user_msg=$(printf '%s' "$user_msg" | jq -Rsa .)

    local payload
    payload=$(cat <<JSONEOF
{
    "messages": [
        {"role": "system", "content": ${system_msg}},
        {"role": "user", "content": ${user_msg}}
    ],
    "temperature": 0.1,
    "max_tokens": 1024,
    "stream": false
}
JSONEOF
)

    local response
    response=$(curl -sf -X POST "$API" \
        -H "Content-Type: application/json" \
        -d "$payload" \
        --max-time 300)

    echo "$response" | jq -r '.choices[0].message.content // "No response from LLM"'
}

SYSTEM_PROMPT="You are a C++23 code reviewer for a resource-constrained shortest path (RCSPP) branch-and-cut solver. The project uses GCC 14, CMake, TBB for parallelism, and optionally HiGHS as MIP solver. Be concise and precise. Focus on correctness, not style."

case "$MODE" in
    build_failure)
        query_llm \
            "$SYSTEM_PROMPT" \
            "The build failed with these compiler errors/warnings. Diagnose the root cause and suggest a minimal fix for each issue. Only suggest changes you are confident about.

Build output:
$CONTEXT"
        ;;

    test_failure)
        query_llm \
            "$SYSTEM_PROMPT" \
            "These tests failed. Analyze the failures, identify likely root causes, and suggest minimal fixes. If you cannot determine the cause, say so.

Test output:
$CONTEXT"
        ;;

    warnings)
        query_llm \
            "$SYSTEM_PROMPT" \
            "The build succeeded but produced these warnings. For each warning, explain whether it could cause bugs and suggest a fix if appropriate. Ignore trivial warnings from third-party code.

Warnings:
$CONTEXT"
        ;;

    review)
        # Collect recent git diff for review
        DIFF=$(cd /workspace && git diff HEAD~1 --stat --diff-filter=ACMR -p -- '*.cpp' '*.h' '*.hpp' | head -500 2>/dev/null || echo "No diff available")

        query_llm \
            "$SYSTEM_PROMPT" \
            "Review the following code changes for bugs, undefined behavior, race conditions, or logic errors. Be concise — only flag real issues, not style nits.

Diff:
$DIFF"
        ;;

    *)
        echo "Unknown mode: $MODE" >&2
        exit 1
        ;;
esac
