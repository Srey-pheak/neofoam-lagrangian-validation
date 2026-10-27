#!/bin/bash
#------------------------------------------------------------------------------
# Test B - forces and drag: runs every NeoFOAM run of test B and collects the CSVs.
#
#   stokes_<scheme>   stokesDrag + gravity, endTime 1 s, output every 0.0005 s
#                     (exact solution for every diameter)
#   sphere_<scheme>   sphereDrag (Putnam, as OpenFOAM) + gravity, endTime 8 s
#   SN_<scheme>       SchillerNaumannDrag + gravity, endTime 8 s
#                     (terminal velocity vs the force balance)
#
# every run with both integrators: <scheme> = RK2 and Euler (6 runs)
#
# Each run is a copy of the reference case in runs/<name>/; its CSV is copied
# to results/<name>.csv.
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/run_B.sh
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the test folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }

# run_one <name> <drag keyword> <endTime> <writeInterval> <scheme>
run_one()
{
    local name=$1 drag=$2 endTime=$3 writeInterval=$4 scheme=$5
    local R=runs/$name
    rm -rf "$R"
    mkdir -p "$R" results
    cp -r 0.orig constant system Allrun Allclean "$R"/

    # the only changes to the reference case:
    sed -i "s/^\(\s*\)sphereDrag;/\1$drag;/" "$R/constant/kinematicCloudProperties"
    sed -i "s/^\(\s*U\s\+\)RK2;/\1$scheme;/" "$R/constant/kinematicCloudProperties"
    sed -i "s/^endTime .*/endTime         $endTime;/" "$R/system/controlDict"
    sed -i "s/^writeInterval .*/writeInterval   $writeInterval;/" "$R/system/controlDict"

    (cd "$R" && ./Allrun > /dev/null)
    if ! tail -3 "$R/log.neoIcoUncoupledKinematicParcelFoam" | grep -q "^End"
    then
        echo "  $name: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"
        return
    fi
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/$name.csv"
    echo "  $name: done -> results/$name.csv"
}

for scheme in RK2 Euler
do
    run_one stokes_$scheme stokesDrag          1 10 $scheme
    run_one sphere_$scheme sphereDrag          8 40 $scheme
    run_one SN_$scheme     SchillerNaumannDrag 8 40 $scheme
done

echo "Plot:  python3 scripts/plot_B.py"
