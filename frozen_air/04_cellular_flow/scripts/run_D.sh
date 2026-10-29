#!/bin/bash
#------------------------------------------------------------------------------
# Test 4 - cellular flow: copies the reference case to runs/D/, runs NeoFOAM
# there, copies the CSV to results/D.csv and the cell values of U (0/U, written by
# setExprFields) to results/U_cells.
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/run_D.sh
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the test folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }

R=runs/D
rm -rf "$R"
mkdir -p "$R" results
cp -r 0.orig constant system Allrun Allclean "$R"/
(cd "$R" && ./Allrun > /dev/null)
if ! tail -3 "$R/log.neoIcoUncoupledKinematicParcelFoam" | grep -q "^End"
then
    echo "FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"
    exit 1
fi
grep "left the domain" "$R/log.neoIcoUncoupledKinematicParcelFoam" | tail -1 || true
cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" results/D.csv
cp "$R/0/U" results/U_cells       # the cell values the solver used (for the same-cell reference)
echo "done -> results/D.csv     Plot:  python3 scripts/plot_D.py | tee results/summary.txt"
