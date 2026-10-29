#!/bin/sh
# Run Test 6: coupled solver, live non-uniform air (lid-driven cavity).
#
# Copies the clean case to runs/<name>/, runs blockMesh + the solver, collects
# the particle CSV into results/<name>.csv and runs the Python replay, which
# checks that the particles follow the air of each step.
#
# The time folders are KEPT in runs/<name>/ (0.005/U ... 1/U): the replay reads
# them. They are gitignored with the rest of runs/.
#
# Usage:
#   scripts/run_F.sh [run-name]        (default: 06_cavity)
#
# Requires OpenFOAM 2406 sourced in the shell and a production build of
# neoIcoKinematicParcelFoam.  Override the binary directory with NEOFOAM_BIN.

set -e

CASE_DIR=$(cd "$(dirname "$0")/.." && pwd)
RUN_NAME=${1:-06_cavity}
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

# ParaView: the particles as lagrangian files in the time folders (as OpenFOAM writes them)
python3 "$CASE_DIR/../write_lagrangian.py" "$RUN_DIR"

# The replay exits non-zero on FAIL: keep going so the summary is still saved.
status=0
python3 "$CASE_DIR/scripts/replay_F.py" "$RUN_DIR" \
    > "$CASE_DIR/results/$RUN_NAME.summary.txt" || status=$?
cat "$CASE_DIR/results/$RUN_NAME.summary.txt"

echo
echo "time folders kept in: $RUN_DIR (needed by the replay, not committed)"
echo "ParaView: open $RUN_DIR/$RUN_NAME.foam"
exit "$status"
