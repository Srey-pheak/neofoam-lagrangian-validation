#!/bin/bash
#------------------------------------------------------------------------------
# Check: bad case input is rejected with a clear error.
#
# Starts from the test A case, breaks ONE thing per run, runs NeoFOAM and
# expects it to stop with a FATAL ERROR. Prints the error message per case.
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/check_input.sh
#------------------------------------------------------------------------------
cd "$(dirname "$0")/.."                          # the checks folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,vader}
NEOFOAM_BIN=${NEOFOAM_BIN:?set NEOFOAM_BIN to <NeoFOAM>/build/production/bin}
CASE=../A_integrators
CP=constant/kinematicCloudProperties

# name | sed command applied to the copy | file
cases=(
  "two_drag_models|s/^\(\s*\)stokesDrag;.*/\1stokesDrag;\n\1sphereDrag;/|$CP"
  "no_drag_model|/^\s*stokesDrag;/d|$CP"
  "unknown_force|s/^\(\s*\)stokesDrag;.*/\1stokesDrag;\n\1pressureGradient;/|$CP"
  "analytical_scheme|s/^\(\s*U\s\+\)Euler;/\1analytical;/|$CP"
  "SOI_not_zero|s/^\(\s*SOI\s\+\)0;/\10.1;/|$CP"
  "injection_outside_mesh|s/^(0.05 0.05 0.05)/(5 5 5)/|constant/kinematicCloudPositions"
)

mkdir -p results
: > results/check_input.txt
pass=0; total=0
for c in "${cases[@]}"
do
    IFS='|' read -r name cmd file <<< "$c"
    R=runs/input_$name
    rm -rf "$R"; mkdir -p "$R"
    cp -r $CASE/0.orig $CASE/constant $CASE/system "$R"/
    sed -i "$cmd" "$R/$file"
    (
        cd "$R"
        restore0Dir > /dev/null
        runApplication blockMesh
        "$NEOFOAM_BIN/neoIcoUncoupledKinematicParcelFoam" > log.neoIcoUncoupledKinematicParcelFoam 2>&1
    )
    total=$((total + 1))
    log=$R/log.neoIcoUncoupledKinematicParcelFoam
    if grep -q "FATAL ERROR" "$log"
    then
        msg=$(grep -A3 "FATAL ERROR" "$log" | grep -v "FATAL ERROR" | grep -v '^\s*$' | head -1 | sed 's/^\s*//')
        pass=$((pass + 1))
        line="  PASS  $name: stopped with: $msg"
    else
        line="  FAIL  $name: did NOT stop with an error (see $log)"
    fi
    echo "$line" | tee -a results/check_input.txt
done
echo "  $pass / $total rejected" | tee -a results/check_input.txt
