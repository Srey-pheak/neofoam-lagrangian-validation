#!/bin/sh
# Run Test 5: coupled solver, uniform carrier flow.
#
# Copies the clean case to runs/<name>/, runs blockMesh + the solver, collects
# the particle CSV into results/<name>.csv and compares it with Test 3.
#
# Usage:
#   scripts/run_E.sh [run-name]        (default: E_uniform)
#
# Requires OpenFOAM 2406 sourced in the shell and a production build of
# neoIcoKinematicParcelFoam.  Override the binary directory with NEOFOAM_BIN.

set -e

CASE_DIR=$(cd "$(dirname "$0")/.." && pwd)
RUN_NAME=${1:-E_uniform}
RUN_DIR="$CASE_DIR/runs/$RUN_NAME"
NEOFOAM_BIN=${NEOFOAM_BIN:?set NEOFOAM_BIN to <NeoFOAM>/build/production/bin}

case "$RUN_NAME" in
    */*|.*)
        echo "invalid run name: $RUN_NAME" >&2
        exit 1
        ;;
esac

case "$RUN_DIR" in
    "$CASE_DIR"/runs/?*)
        ;;
    *)
        echo "refusing to delete unexpected path: $RUN_DIR" >&2
        exit 1
        ;;
esac

if ! command -v blockMesh >/dev/null 2>&1; then
    echo "blockMesh not found. Source OpenFOAM first, e.g.:" >&2
    echo "  source /usr/lib/openfoam/openfoam2406/etc/bashrc" >&2
    exit 1
fi

if [ ! -x "$NEOFOAM_BIN/neoIcoKinematicParcelFoam" ]; then
    echo "solver not found: $NEOFOAM_BIN/neoIcoKinematicParcelFoam" >&2
    echo "build it or set NEOFOAM_BIN" >&2
    exit 1
fi

echo "case:   $CASE_DIR"
echo "run:    $RUN_DIR"
echo "solver: $NEOFOAM_BIN/neoIcoKinematicParcelFoam"

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cp -r "$CASE_DIR"/0.orig "$CASE_DIR"/constant "$CASE_DIR"/system \
      "$CASE_DIR"/Allrun "$CASE_DIR"/Allclean "$RUN_DIR"/

cd "$RUN_DIR"
NEOFOAM_BIN="$NEOFOAM_BIN" sh ./Allrun

mkdir -p "$CASE_DIR/results"
cp "$RUN_DIR"/postProcessing/lagrangian/kinematicCloud/particles.csv \
   "$CASE_DIR/results/$RUN_NAME.csv"

python3 "$CASE_DIR/scripts/compare_to_C.py" "$CASE_DIR/results/$RUN_NAME.csv" \
    > "$CASE_DIR/results/$RUN_NAME.summary.txt"
cat "$CASE_DIR/results/$RUN_NAME.summary.txt"
