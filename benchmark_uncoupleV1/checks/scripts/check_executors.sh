#!/bin/bash
#------------------------------------------------------------------------------
# Check: CPU and GPU executors give the same result as Serial.
#
# Runs the test C case (3 particles, 1000 cells, crossing cells, one leaving)
# with -executor Serial, CPU and (if available) GPU, and compares the CSVs.
# Expected: bitwise identical for Serial vs CPU; GPU identical or round-off.
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/check_executors.sh [GPU]
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the checks folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,vader}
NEOFOAM_BIN=${NEOFOAM_BIN:?set NEOFOAM_BIN to <NeoFOAM>/build/production/bin}
CASE=../C_3dBoundaries

executors="Serial CPU"
[ "$1" = GPU ] && executors="$executors GPU"

mkdir -p results
for exe in $executors
do
    R=runs/executor_$exe
    rm -rf "$R"; mkdir -p "$R"
    cp -r $CASE/0.orig $CASE/constant $CASE/system "$R"/
    (
        cd "$R"
        restore0Dir > /dev/null
        runApplication blockMesh
        "$NEOFOAM_BIN/neoIcoUncoupledKinematicParcelFoam" -executor $exe > log.neoIcoUncoupledKinematicParcelFoam 2>&1
        foamListTimes -rm > /dev/null
    ) || { echo "  $exe: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"; continue; }
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/executor_$exe.csv"
    echo "  $exe: done -> results/executor_$exe.csv"
done

echo
for exe in $executors
do
    [ "$exe" = Serial ] && continue
    [ -f "results/executor_$exe.csv" ] || continue
    if cmp -s results/executor_Serial.csv "results/executor_$exe.csv"
    then
        echo "  Serial vs $exe: bitwise identical  ->  PASS"
    else
        python3 - "$exe" <<'PY'
import sys, numpy as np
a = np.genfromtxt("results/executor_Serial.csv", delimiter=",", skip_header=1)
b = np.genfromtxt(f"results/executor_{sys.argv[1]}.csv", delimiter=",", skip_header=1)
if a.shape != b.shape:
    print(f"  Serial vs {sys.argv[1]}: different number of rows  ->  FAIL")
else:
    d = np.abs(a - b).max()
    print(f"  Serial vs {sys.argv[1]}: max difference {d:.2e}  ->  {'PASS (round-off)' if d < 1e-12 else 'FAIL'}")
PY
    fi
done
