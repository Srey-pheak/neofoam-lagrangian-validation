#!/bin/bash
#------------------------------------------------------------------------------
# Reference runs: the Benchmark 1 cases with OpenFOAM v2406
# icoUncoupledKinematicParcelFoam, and the IDENTICAL case with NeoFOAM, for a
# one-to-one comparison.
#
# OpenFOAM has only sphereDrag (Putnam) and Euler of our models (no stokesDrag,
# no SchillerNaumann, no RK2), so every reference case uses sphereDrag + Euler:
#
#   A_<dt>   test 1 case (box 0.1 m, 1 cell, 1 particle), dt = 2e-3 ... 1.25e-4
#   B        test 2 case (35 m column, 1 cell, 5 sizes), dt = 5e-5, 8 s
#   C        test 3 case (1 m box, 10x10x10 cells, wind, 3 particles), dt = 5e-4, 0.3 s
#
# The two codes run separately, each in its own run folder:
#   runs/openfoam/<name>/  ->  results/openfoam/<name>.csv  (read_openfoam.py)
#   runs/neofoam/<name>/   ->  results/neofoam/<name>.csv
#
# Usage (OpenFOAM loaded, from anywhere):
#   scripts/run_reference.sh            both codes, all cases
#   scripts/run_reference.sh openfoam   OpenFOAM only
#   scripts/run_reference.sh neofoam    NeoFOAM only
#   scripts/run_reference.sh both A_dt1e-3   one case only (name prefix: A, A_dt1e-3, B, C)
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the reference_openfoam folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,vader}
NEOFOAM_BIN=${NEOFOAM_BIN:?set NEOFOAM_BIN to <NeoFOAM>/build/production/bin}
WHICH=${1:-both}
ONLY=${2:-}

# make_case <run folder> <test folder> <dt> <endTime> <writeInterval> <application>
make_case()
{
    local R=$1 T=../$2
    rm -rf "$R"; mkdir -p "$R"
    cp -r "$T/0.orig" "$T/constant" "$T/system" "$R"/
    # the only changes to the test case: sphereDrag + Euler, time step, end time, output
    sed -i -E "s/^(\s*)(stokesDrag|sphereDrag|SchillerNaumannDrag);.*/\1sphereDrag;/" "$R/constant/kinematicCloudProperties"
    sed -i -E "s/^(\s*U\s+)(Euler|RK2);/\1Euler;/" "$R/constant/kinematicCloudProperties"
    sed -i "s/^deltaT .*/deltaT          $3;/; s/^endTime .*/endTime         $4;/; s/^writeInterval .*/writeInterval   $5;/" "$R/system/controlDict"
    sed -i "s/^application .*/application     $6;/" "$R/system/controlDict"
}

# run_openfoam <name> <test folder> <dt> <endTime> <writeInterval>
run_openfoam()
{
    local R=runs/openfoam/$1
    make_case "$R" "$2" "$3" "$4" "$5" icoUncoupledKinematicParcelFoam
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh && runApplication icoUncoupledKinematicParcelFoam) \
        || { echo "  openfoam $1: FAILED, see $R/log.icoUncoupledKinematicParcelFoam"; return; }
    mkdir -p results/openfoam
    python3 scripts/read_openfoam.py "$R" "results/openfoam/$1.csv" | sed "s/^/  openfoam $1: /"
}

# run_neofoam <name> <test folder> <dt> <endTime> <writeInterval>
run_neofoam()
{
    local R=runs/neofoam/$1
    make_case "$R" "$2" "$3" "$4" "$5" neoIcoUncoupledKinematicParcelFoam
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh \
        && runApplication "$NEOFOAM_BIN/neoIcoUncoupledKinematicParcelFoam" && foamListTimes -rm > /dev/null) \
        || { echo "  neofoam $1: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"; return; }
    mkdir -p results/neofoam
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/neofoam/$1.csv"
    echo "  neofoam $1: done -> results/neofoam/$1.csv"
}

run_both()
{
    [ -z "$ONLY" ] || [[ "$1" == "$ONLY"* ]] || return 0
    [ "$WHICH" = neofoam ] || run_openfoam "$@"
    [ "$WHICH" = openfoam ] || run_neofoam "$@"
}

echo "A: integrators case, sphereDrag + Euler, dt sweep"
for dt in 2e-3 1e-3 5e-4 2.5e-4 1.25e-4
do
    run_both "A_dt$dt" 01_integrators $dt 0.2 1
done

echo "B: settling column, sphereDrag + Euler"
run_both B 02_drag_settling 5e-5 8 40

echo "C: 3D box with wind, sphereDrag + Euler"
run_both C 03_motion_3d_exit 5e-4 0.3 1

echo "Compare:  python3 scripts/compare.py"
