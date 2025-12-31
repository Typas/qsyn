#! /usr/bin/env bash
set -euo pipefail

# Ensure clean build directory to avoid CMake cache issues
# CMake caches compiler choice in CMakeCache.txt, so we need to remove
# the build directory to ensure we use the compiler specified in ENV
rm -rf /app/build

# Explicitly specify compiler to ensure deterministic builds
# This overrides any CMake cache and ensures we use the compiler
# specified in the Dockerfile (gcc/g++ or clang/clang++)
cmake -B /app/build -S /app/qsyn \
  -DCMAKE_C_COMPILER="${CC:-/usr/bin/gcc}" \
  -DCMAKE_CXX_COMPILER="${CXX:-/usr/bin/g++}"

# Verify compiler selection (for debugging)
echo "Using compiler:"
grep -E "CMAKE_(C|CXX)_COMPILER:" /app/build/CMakeCache.txt || true

cmake --build /app/build --parallel "$(nproc)"

cd /app/qsyn || exit 1
/app/build/qsyn-unit-test || exit 1
/app/qsyn/scripts/RUN_TESTS --qsyn /app/build/qsyn "$@"
