#!/bin/bash
set -e # Stop script on first error

BUILD_DIR="build"

echo "--- Configuring and building project ---"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake -D compiler=intel ..
make -j$(nproc) # Use multiple cores for faster build

echo "--- Running solver ---"
./LaplaceSolver