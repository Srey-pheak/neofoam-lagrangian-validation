#!/bin/sh
# Run Test 6 on two executors and compare the particle CSVs.
#
# Default: Serial (CPU build) vs GPU (CUDA build). Each run uses its own solver
# binary, because a NeoFOAM build has one backend:
#
#   $NEOFOAM_BUILD/production/bin       Kokkos Serial  (no OpenMP, no CUDA)
#   $NEOFOAM_BUILD/production-cuda/bin  Kokkos Serial + CUDA
#
# Both runs are checked by scripts/replay_F.py; then the CSVs are compared.
# Expected: identical to round-off (the GPUs and the CPU do the same float64
# arithmetic in a different order). A large difference is a real bug, not
# "GPU noise".
#
# Usage:
#   scripts/run_F_executors.sh            # Serial vs GPU
#   scripts/run_F_executors.sh Serial     # just Serial
#
# Requires OpenFOAM sourced and NEOFOAM_BUILD=<NeoFOAM>/build.

set -e

CASE_DIR=$(cd "$(dirname "$0")/.." && pwd)
NEOFOAM_BUILD=${NEOFOAM_BUILD:?set NEOFOAM_BUILD to <NeoFOAM>/build}

list=${1:-"Serial GPU"}

if ! command -v blockMesh >/dev/null 2>&1; then
    echo "blockMesh not found. Source OpenFOAM first, e.g.:" >&2
    echo "  source /usr/lib/openfoam/openfoam2406/etc/bashrc" >&2
    exit 1
fi

for exe in $list
do
    case "$exe" in
        Serial) bin="$NEOFOAM_BUILD/production/bin" ;;
        GPU)    bin="$NEOFOAM_BUILD/production-cuda/bin" ;;
        *) echo "unknown executor: $exe (use Serial or GPU)" >&2; exit 1 ;;
    esac

    if [ ! -x "$bin/neoIcoKinematicParcelFoam" ]; then
        echo "missing $bin/neoIcoKinematicParcelFoam" >&2
        if [ "$exe" = GPU ]; then
            echo "  the CUDA build must include the coupled solver; rebuild it:" >&2
            echo "    source /usr/lib/openfoam/openfoam2406/etc/bashrc" >&2
            echo "    export PATH=/usr/local/cuda/bin:\$PATH" >&2
            echo "    cd <NeoFOAM> && NEON_DEVICE=-cuda cmake --preset production" >&2
            echo "    NEON_DEVICE=-cuda cmake --build --preset production -j 8" >&2
        fi
        exit 1
    fi

    echo "=== $exe: $bin/neoIcoKinematicParcelFoam ==="
    start=$(date +%s)
    NEOFOAM_BIN="$bin" NEOFOAM_EXECUTOR="$exe" \
        "$CASE_DIR/scripts/run_F.sh" "F_$exe" | tail -20
    end=$(date +%s)
    echo "  wall time: $((end - start)) s"
    echo
done

# ─── compare the executor runs with each other ───
[ -f "$CASE_DIR/results/F_Serial.csv" ] || exit 0
for exe in $list
do
    [ "$exe" = Serial ] && continue
    other="$CASE_DIR/results/F_$exe.csv"
    [ -f "$other" ] || continue
    if cmp -s "$CASE_DIR/results/F_Serial.csv" "$other"
    then
        echo "Serial vs $exe: bitwise identical  ->  PASS"
    else
        python3 - "$CASE_DIR/results/F_Serial.csv" "$other" "$exe" <<'PY'
import sys
import numpy as np
a = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
b = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
if a.dtype.names != b.dtype.names or len(a) != len(b):
    print(f"Serial vs {sys.argv[3]}: different shape  ->  FAIL")
    sys.exit(1)
worst = max(float(np.abs(a[c] - b[c]).max()) for c in a.dtype.names)
print(f"Serial vs {sys.argv[3]}: max |difference| = {worst:.2e} over all columns  "
      f"->  {'PASS (round-off)' if worst < 1e-12 else 'FAIL'}")
PY
    fi
done
