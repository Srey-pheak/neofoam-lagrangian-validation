#!/bin/bash
#------------------------------------------------------------------------------
# Reference runs for test 6 (lid-driven cavity, air solved every step): the
# same case with OpenFOAM v2406 and with NeoFOAM, for a one-to-one comparison.
#
# Changes to the test 6 case, the same in every run (make_case):
#   - sphereDrag + Euler: OpenFOAM has no stokesDrag and no RK2 (as the
#     Benchmark 1 reference)
#   - linear solver tolerance 1e-12: NeoN and OpenFOAM normalise residuals
#     differently, so at the test 6 tolerances (1e-6, 1e-5) the two airs differ
#     by the solver convergence (~2e-4 m/s), not by the method
#
# Runs:
#   neofoam/F          NeoFOAM neoIcoKinematicParcelFoam: air + particles
#   openfoam/F_fluid   OpenFOAM icoFoam: the same air equations as neoIcoFoam
#   openfoam/F_chain   OpenFOAM particles (icoUncoupledKinematicParcelFoam) in
#                      the NeoFOAM air, one time step per launch: before the step
#                      from t_k the NeoFOAM U(t_k) is copied into t_k/U. The first
#                      step runs in the start air U(0) = 0, in which OpenFOAM's
#                      double move of a just-injected particle changes nothing; so
#                      OpenFOAM at t_k+1 = NeoFOAM at t_k, moved in the same air.
#   openfoam/F         information only: OpenFOAM's own one-way solver
#                      kinematicParcelFoam (pimpleFoam air, particles moved BEFORE
#                      the air is solved)
#
#   runs/<code>/<name>/  ->  results/<code>/<name>.csv
#
# ParaView: every run folder gets a <code>_<name>.foam file (open it with the
# built-in OpenFOAM reader); the NeoFOAM particles are also written as
# <time>/lagrangian/kinematicCloud/ (../write_lagrangian.py), as OpenFOAM writes them.
#
# Usage (OpenFOAM loaded, from anywhere):
#   scripts/run_reference.sh            both codes, all runs
#   scripts/run_reference.sh openfoam   OpenFOAM only (needs runs/neofoam/F)
#   scripts/run_reference.sh neofoam    NeoFOAM only
#------------------------------------------------------------------------------
set -e
cd "$(dirname "$0")/.."                          # the reference_openfoam folder
[ -n "$WM_PROJECT_VERSION" ] || { echo "Load OpenFOAM first: source /usr/lib/openfoam/openfoam2406/etc/bashrc"; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,vader}
WHICH=${1:-both}
T=../06_cavity                                    # the test 6 case

# make_case <run folder> <application>
make_case()
{
    local R=$1
    rm -rf "$R"; mkdir -p "$R"
    cp -r "$T/0.orig" "$T/constant" "$T/system" "$R"/
    sed -i -E "s/^(\s*)(stokesDrag|sphereDrag|SchillerNaumannDrag);.*/\1sphereDrag;/" "$R/constant/kinematicCloudProperties"
    sed -i -E "s/^(\s*U\s+)(Euler|RK2);/\1Euler;/" "$R/constant/kinematicCloudProperties"
    sed -i -E "s/^(\s*tolerance\s+).*/\11e-12;/" "$R/system/fvSolution"
    sed -i "s/^application .*/application     $2;/" "$R/system/controlDict"
    # OpenFOAM's PISO/PIMPLE solve the last corrector with the <field>Final solver
    # (regex keys: OpenFOAM only, NeoFOAM reads the plain p and U entries)
    if [ "$2" != neoIcoKinematicParcelFoam ]; then
        sed -i -E 's/^(\s*)p$/\1"(p|pFinal)"/; s/^(\s*)U$/\1"(U|UFinal)"/' "$R/system/fvSolution"
    fi
}

# OpenFOAM kinematicParcelFoam: PIMPLE (as PISO), laminar stress schemes, and the
# name of the carrier viscosity field (muc) for the cloud interpolation
add_pimple()
{
    local R=$1
    cat >> "$R/system/fvSolution" <<'EOF'

PIMPLE
{
    momentumPredictor        yes;
    nOuterCorrectors         1;
    nCorrectors              2;
    nNonOrthogonalCorrectors 0;
    pRefCell                 0;
    pRefValue                0;
}
EOF
    sed -i -E 's/^(\s*)(div\(phi,U\).*)$/\1\2\n\1div((nuEff*dev2(T(grad(U))))) Gauss linear;/' "$R/system/fvSchemes"
    sed -i -E 's/^(\s*)(laplacian\(nu,U\).*)$/\1\2\n\1laplacian(nuEff,U)      Gauss linear uncorrected;/' "$R/system/fvSchemes"
    sed -i -E 's/^(\s*)mu(\s+)cell;/\1mu\2cell;\n\1muc\2cell;/' "$R/constant/kinematicCloudProperties"
}

run_neofoam()
{
    NEOFOAM_BIN=${NEOFOAM_BIN:?set NEOFOAM_BIN to <NeoFOAM>/build/production/bin}
    local R=runs/neofoam/F
    make_case "$R" neoIcoKinematicParcelFoam
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh \
        && runApplication "$NEOFOAM_BIN/neoIcoKinematicParcelFoam") \
        || { echo "  neofoam F: FAILED, see $R/log.neoIcoKinematicParcelFoam"; return; }
    mkdir -p results/neofoam
    cp "$R/postProcessing/lagrangian/kinematicCloud/particles.csv" "results/neofoam/F.csv"
    python3 ../write_lagrangian.py "$R" > /dev/null          # particles as lagrangian files, for ParaView
    touch "$R/neofoam_F.foam"
    echo "  neofoam F: done -> results/neofoam/F.csv (U every step kept in $R)"
}

run_openfoam_fluid()
{
    local R=runs/openfoam/F_fluid
    make_case "$R" icoFoam
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh && runApplication icoFoam) \
        || { echo "  openfoam F_fluid: FAILED, see $R/log.icoFoam"; return; }
    touch "$R/openfoam_F_fluid.foam"
    echo "  openfoam F_fluid: done -> $R (U every step)"
}

run_openfoam_chain()
{
    local N=runs/neofoam/F R=runs/openfoam/F_chain
    [ -f "$N/3/U" ] || { echo "  openfoam F_chain: needs the NeoFOAM run first ($N)"; return; }
    make_case "$R" icoUncoupledKinematicParcelFoam
    sed -i "s/^startFrom .*/startFrom       latestTime;/" "$R/system/controlDict"
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh) \
        || { echo "  openfoam F_chain: blockMesh FAILED"; return; }
    local times=($(cd "$N" && foamListTimes -withZero))
    local k
    for ((k = 0; k < ${#times[@]} - 1; k++))
    do
        cp "$N/${times[k]}/U" "$R/${times[k]}/U"            # the air NeoFOAM used for this step
        sed -i "s/^endTime .*/endTime         ${times[k+1]};/" "$R/system/controlDict"
        (cd "$R" && icoUncoupledKinematicParcelFoam >> log.icoUncoupledKinematicParcelFoam 2>&1) \
            || { echo "  openfoam F_chain: FAILED at t = ${times[k]}, see $R/log.icoUncoupledKinematicParcelFoam"; return; }
    done
    mkdir -p results/openfoam
    touch "$R/openfoam_F_chain.foam"
    python3 scripts/read_openfoam.py "$R" "results/openfoam/F_chain.csv" | sed "s/^/  openfoam F_chain: /"
}

run_openfoam_parcel()
{
    local R=runs/openfoam/F
    make_case "$R" kinematicParcelFoam
    add_pimple "$R"
    (cd "$R" && restore0Dir > /dev/null && runApplication blockMesh && runApplication kinematicParcelFoam) \
        || { echo "  openfoam F: FAILED, see $R/log.kinematicParcelFoam"; return; }
    mkdir -p results/openfoam
    touch "$R/openfoam_F.foam"
    python3 scripts/read_openfoam.py "$R" "results/openfoam/F.csv" | sed "s/^/  openfoam F: /"
}

echo "F: lid-driven cavity, air solved every step, sphereDrag + Euler, tolerance 1e-12"
[ "$WHICH" = openfoam ] || run_neofoam
if [ "$WHICH" != neofoam ]; then
    run_openfoam_fluid
    run_openfoam_chain
    run_openfoam_parcel
fi

echo "Compare:  python3 scripts/compare.py | tee results/summary.txt"
