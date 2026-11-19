#!/bin/bash
# Clone and build PathWyse into benchmarks/.pathwyse/, pinned to one release.
#
# The pin is the 1.1.0 tag and the commit it names; a checkout at any other
# commit is refused rather than reused, so results always come from the
# pinned revision. To move the pin, change both lines and delete .pathwyse/.
#
# The paper's PathWyse rows (benchmarks/pathwyse.csv, dated 2026-03-07) came
# from an unpinned clone of upstream main. Main was at 0.2.0 (d53c01b) from
# 2025-06-13 to 2026-07-09, so that is most likely the revision behind them;
# an older clone (0.1.4, 397b871) cannot be ruled out, since this script used
# to reuse any existing checkout.
set -euo pipefail

PATHWYSE_TAG="1.1.0"
PATHWYSE_COMMIT="15a7468fc13891d77bf25bbc9235a07101f271ba"

DIR="$(cd "$(dirname "$0")" && pwd)"
PW_DIR="$DIR/.pathwyse"

if [ -d "$PW_DIR" ]; then
    have="$(git -C "$PW_DIR" rev-parse HEAD)"
    if [ "$have" != "$PATHWYSE_COMMIT" ]; then
        echo "Error: $PW_DIR is at $have, not PathWyse $PATHWYSE_TAG ($PATHWYSE_COMMIT)." >&2
        echo "Delete it and run this script again." >&2
        exit 1
    fi
    if [ -n "$(git -C "$PW_DIR" status --porcelain --untracked-files=no)" ]; then
        echo "Error: $PW_DIR has local changes; delete it and run this script again." >&2
        exit 1
    fi
    # Build anyway: incremental, so a no-op for a binary built at this
    # commit, and it rebuilds one left from a checkout moved to the tag.
    echo "PathWyse $PATHWYSE_TAG already cloned at $PW_DIR"
else
    echo "Cloning PathWyse $PATHWYSE_TAG..."
    git -c advice.detachedHead=false clone --depth 1 --branch "$PATHWYSE_TAG" \
        https://github.com/pathwyse/pathwyse.git "$PW_DIR"
    have="$(git -C "$PW_DIR" rev-parse HEAD)"
    if [ "$have" != "$PATHWYSE_COMMIT" ]; then
        echo "Error: tag $PATHWYSE_TAG resolves to $have, expected $PATHWYSE_COMMIT." >&2
        exit 1
    fi
fi

echo "Building PathWyse..."
cd "$PW_DIR"
cmake -B build -DCMAKE_BUILD_TYPE=Release
if command -v nproc >/dev/null; then
    NCPU="$(nproc)"
else
    NCPU="$(sysctl -n hw.logicalcpu)"
fi
cmake --build build -j"$NCPU"

echo "PathWyse $PATHWYSE_TAG built successfully at $PW_DIR/bin/pathwyse"
