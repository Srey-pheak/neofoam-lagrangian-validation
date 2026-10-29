#!/bin/bash
#------------------------------------------------------------------------------
# CPU vs GPU: runs the reference cases of tests 1, 2 and 3 (unchanged) with
#   CPU   build/production       -executor Serial   (as all other tests)
#   GPU   build/production-cuda  -executor GPU
# and copies each CSV to results/<case>_<CPU|GPU>.csv. The run time of each
# run is written to results/runtime.txt.
#
#   A   01_integrators reference case (Euler, dt = 1e-3, 1 particle)
#   B   02_drag_settling  reference case (sphereDrag + RK2, 5 particles, 160 000 steps)
#   C   03_motion_3d_exit reference case (stokesDrag + RK2, 3 particles, 1000 cells)
#
# Usage (OpenFOAM loaded, from anywhere):   scripts/run_cpu_gpu.sh [A|B|C]
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the cpu_vs_gpu folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,vader}
BUILD=${NEOFOAM_BUILD:?set NEOFOAM_BUILD to <NeoFOAM>/build (needs production/ and production-cuda/)}
CASES=${1:-A B C}

# run_one <case> <test folder> <CPU|GPU>
run_one()
{
    local name=$1 T=../$2 dev=$3 bin exe
    if [ "$dev" = GPU ]; then bin=$BUILD/production-cuda/bin; exe=GPU; else bin=$BUILD/production/bin; exe=Serial; fi
    local R=runs/${name}_$dev
    rm -rf "$R"; mkdir -p "$R" results
    cp -r "$T/0.orig" "$T/constant" "$T/system" "$R"/
    local start=$(date +%s.%N)
    (
        cd "$R"
        restore0Dir > /dev/null
        runApplication blockMesh
        "$bin/neoIcoUncoupledKinematicParcelFoam" -executor $exe > log.neoIcoUncoupledKinematicParcelFoam 2>&1
        foamListTimes -rm > /dev/null
    ) || { echo "  ${name}_$dev: FAILED, see $R/log.neoIcoUncoupledKinematicParcelFoam"; return; }
    local sec=$(echo "$(date +%s.%N) - $start" | bc)
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/${name}_$dev.csv"
    sed -i "/^${name}_$dev /d" results/runtime.txt 2>/dev/null || true
    printf "%s %.1f\n" "${name}_$dev" "$sec" >> results/runtime.txt
    printf "  %s: done in %.1f s -> results/%s.csv\n" "${name}_$dev" "$sec" "${name}_$dev"
}

for c in $CASES
do
    case $c in
        A) T=01_integrators ;;
        B) T=02_drag_settling ;;
        C) T=03_motion_3d_exit ;;
    esac
    run_one $c $T CPU
    run_one $c $T GPU
done

echo "Plot:  python3 scripts/plot_cpu_gpu.py | tee results/summary.txt"
