#!/bin/sh
# Compile the GUIDANCE3 helper programs from source into the directory given as $1.
# Used by the Dockerfile for platforms without a bundled binary (linux/arm64);
# the bundled x86-64 Linux and macOS binaries in guidance3/programs/ are built from the same sources.
set -e
OUT=$(realpath "$1")
SRC=$(dirname "$(realpath "$0")")
JOBS=$(nproc 2>/dev/null || echo 4)
mkdir -p "$OUT"

# GUIDANCE2 tools: msa_set_score, isEqualTree, removeTaxa
make -C "$SRC/guidance2/libs" -j"$JOBS"
for p in isEqualTree removeTaxa msa_set_score; do
    make -C "$SRC/guidance2/programs/$p" -j"$JOBS"
    cp "$SRC/guidance2/programs/$p/$p" "$OUT/"
done

# features_for_msas (alignment features for the DL model)
cmake -S "$SRC/features_for_msas" -B /tmp/features_build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/features_build -j"$JOBS"
cp /tmp/features_build/features_for_msas "$OUT/"
