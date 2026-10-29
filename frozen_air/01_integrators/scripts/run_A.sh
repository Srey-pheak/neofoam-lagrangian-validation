#!/bin/bash
#------------------------------------------------------------------------------
# Test 1 - integrators: runs every NeoFOAM run of test 1 and collects the CSVs.
#
#   1. order of accuracy: Euler and RK2, dt = 2e-3 ... 1.25e-4 (each half the previous)
#      on the reference case (box 0.1 m, endTime 0.2 s)
#   2. size of the time step: Euler and RK2 at a small, medium and big dt
#      (0.001, 0.02, 0.05 s = 0.05, 1, 2.6 tau_p) in a big box (10 m, particle at the
#      centre), so a run-away RK2 particle stays inside
#
# Each run is a copy of the reference case in runs/<name>/ (the reference case
# itself is never run). The CSV of each run is copied to
#   results/<scheme>_dt<dt>.csv            (part 1)
#   results/timestep/<scheme>_dt<dt>.csv   (part 2)
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/run_A.sh
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the test folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }

# run_one <name> <scheme> <dt> <results-file> [big]
run_one()
{
    local name=$1 scheme=$2 dt=$3 out=$4 big=$5
    local R=runs/$name
    rm -rf "$R"
    mkdir -p "$R" "$(dirname "$out")"
    cp -r 0.orig constant system Allrun Allclean "$R"/

    # the only changes to the reference case:
    sed -i "s/^deltaT .*/deltaT          $dt;/" "$R/system/controlDict"
    sed -i "s/^\(\s*U\s\+\)Euler;/\1$scheme;/" "$R/constant/kinematicCloudProperties"
    if [ "$big" = big ]
    then
        sed -i "s/^L .*/L       10;/" "$R/system/blockMeshDict"
        sed -i "s/^(0.05 0.05 0.05)/(5 5 5)/" "$R/constant/kinematicCloudPositions"
    fi

    (cd "$R" && ./Allrun > /dev/null)
    if ! tail -3 "$R/log.neoIcoUncoupledKinematicParcelFoam" | grep -q "^End"
    then
        echo "  $name: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"
        return
    fi
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "$out"
    echo "  $name: done -> $out"
}

echo "1. Order of accuracy"
for scheme in Euler RK2
do
    for dt in 2e-3 1e-3 5e-4 2.5e-4 1.25e-4
    do
        run_one "${scheme}_dt$dt" $scheme $dt "results/${scheme}_dt$dt.csv"
    done
done

echo "2. Size of the time step (big box)"
for scheme in Euler RK2
do
    for dt in 0.001 0.02 0.05
    do
        run_one "timestep_${scheme}_dt$dt" $scheme $dt "results/timestep/${scheme}_dt$dt.csv" big
    done
done

echo "Plot:  python3 scripts/plot_A.py"
