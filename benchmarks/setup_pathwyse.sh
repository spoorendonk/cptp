#!/bin/bash
# Clone and build PathWyse into benchmarks/.pathwyse/
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
PW_DIR="$DIR/.pathwyse"

if [ -d "$PW_DIR" ]; then
    echo "PathWyse already set up at $PW_DIR"
    exit 0
fi

echo "Cloning PathWyse..."
git clone https://github.com/pathwyse/pathwyse.git "$PW_DIR"

echo "Building PathWyse..."
cd "$PW_DIR"
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

echo "PathWyse built successfully at $PW_DIR/bin/pathwyse"
