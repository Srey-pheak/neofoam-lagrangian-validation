#!/bin/bash
#------------------------------------------------------------------------------
# Test C - 3D motion and boundaries: for each integrator (RK2, Euler) copies the
# reference case to runs/C_<scheme>/, runs NeoFOAM there, copies the CSV to
# results/C_<scheme>.csv (2 runs).
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/run_C.sh
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the test folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }

for scheme in RK2 Euler
do
    R=runs/C_$scheme
    rm -rf "$R"
    mkdir -p "$R" results
    cp -r 0.orig constant system Allrun Allclean "$R"/

    # the only change to the reference case:
    sed -i "s/^\(\s*U\s\+\)RK2;/\1$scheme;/" "$R/constant/kinematicCloudProperties"

    (cd "$R" && ./Allrun > /dev/null)
    if ! tail -3 "$R/log.neoIcoUncoupledKinematicParcelFoam" | grep -q "^End"
    then
        echo "  C_$scheme: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"
        exit 1
    fi
    grep "left the domain" "$R/log.neoIcoUncoupledKinematicParcelFoam" | head -1 || true
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/C_$scheme.csv"
    echo "  C_$scheme: done -> results/C_$scheme.csv"
done
echo "Plot:  python3 scripts/plot_C.py"
