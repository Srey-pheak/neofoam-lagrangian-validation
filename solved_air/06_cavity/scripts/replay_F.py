#!/usr/bin/env python3
"""Test 6 — lid-driven cavity, air solved every step (neoIcoKinematicParcelFoam).

The air starts at rest and the moving lid drags it into a swirl, so every cell
value of U changes with time and with position. The particles are moved by the
air of the steps, one-way. This script repeats the particle steps in Python,

    step n:  cell = the cell of the particle  (arithmetic, uniform 20 x 20 mesh)
             uf   = U(t_n) in that cell       <- the air NeoFOAM wrote at THIS step
             one RK2 step, Stokes drag, no gravity

and compares with NeoFOAM's CSV. If the solver hands the new air to the
particles in every step, the two agree to round-off.

Two controls make the test able to fail: the same replay with the wrong air.
  control 1: U(t_n-1), one step too old   -> must differ by more than 1e-6 m
  control 2: U(0) = 0, the start air      -> must differ by more than 1e-6 m

Pass: max |x_NeoFOAM - x_Python| < 1e-9 m for the correct air, and both
      controls above 1e-6 m.

Usage:
    python3 scripts/replay_F.py [runs/F]

Defaults to runs/F, i.e. the folder made by scripts/run_F.sh.
"""

import sys
from pathlib import Path

import numpy as np

# ─── Case parameters: must match the case files ───
NX = NY = 20                # system/blockMeshDict
LX = LY = 1.0               # m
LZ = 0.1                    # m, single cell (frontAndBack is empty)
RHO_P = 2500.0              # constantProperties/rho0
RHO_F = 1.2                 # transportProperties/rhoInf
NU = 0.01                   # transportProperties/nu
MU = RHO_F * NU             # = 0.012 Pa s
DT = 5e-3                   # system/controlDict/deltaT
D_P = 1e-3                  # injectionModels/particle_d1mm

DX = LX / NX
TOL_X = 1e-9                # m, correct air
TOL_CONTROL = 1e-6          # m, the wrong air must be at least this far off

# Stokes: F = (pi/8) mu d (CdRe) slip with CdRe = 24 (NeoFOAM DragForce.hpp),
# m = rho_p pi/6 d^3  ->  k/m = 18 mu / (rho_p d^2) = 1/tau_p.
K_OVER_M = np.pi / 8.0 * MU * D_P * 24.0 / (RHO_P * np.pi / 6.0 * D_P**3)


def read_U(path):
    """Cell values of U from one time folder, shape (NX*NY, 3)."""
    text = Path(path).read_text()
    body = text.split("internalField", 1)[1].split("boundaryField")[0]

    if body.lstrip().startswith("uniform"):
        v = body[body.index("(") + 1 : body.index(")")].split()
        return np.tile(np.array(v, float), (NX * NY, 1))

    # nonuniform List<vector>: "<count>" then one "(x y z)" per line
    rows = []
    for line in body.splitlines():
        s = line.strip()
        if s.startswith("(") and s.endswith(")"):
            vals = s[1:-1].split()
            if len(vals) == 3:                      # skip the lone list bracket
                rows.append(vals)
    if not rows:
        raise SystemExit(f"no cell values found in {path}")
    return np.array(rows, float)


def load_particles(csv_path):
    """NeoFOAM CSV -> (X [time, particle, 3], times [time], d [particle])."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n = len(np.unique(data["id"]))
    n_t = len(data) // n
    X = np.stack([data[c][: n * n_t].reshape(n_t, n) for c in "xyz"], axis=2)
    t = data["t"][: n * n_t].reshape(n_t, n)[:, 0]
    return X, t, data["d"][:n]


def inside(x):
    return (
        (x[:, 0] >= 0) & (x[:, 0] < LX)
        & (x[:, 1] >= 0) & (x[:, 1] < LY)
        & (x[:, 2] >= 0) & (x[:, 2] < LZ)
    )


def replay(U, X):
    """Air array 'U[step]' -> particle positions [time, particle, 3]."""
    x = X[0].copy()
    u = np.zeros_like(x)                            # released at rest: U0 = (0 0 0)
    active = np.ones(len(x), bool)
    out = [x.copy()]

    for step in range(1, min(len(U), len(X))):
        active &= inside(x)

        # uniform mesh: cellID = i + NX*j (i fastest, single cell in z)
        i = np.clip(np.floor(x[:, 0] / DX).astype(int), 0, NX - 1)
        j = np.clip(np.floor(x[:, 1] / DX).astype(int), 0, NY - 1)
        uf = U[step][i + NX * j]

        a1 = (uf - u) * K_OVER_M
        u_star = u + a1 * DT
        a2 = (uf - u_star) * K_OVER_M              # same uf: one locate per step
        u_new = u + (a1 + a2) * (0.5 * DT)
        x_new = x + (u + u_star) * (0.5 * DT)

        u = np.where(active[:, None], u_new, u)
        x = np.where(active[:, None], x_new, x)
        out.append(x.copy())

    return np.array(out)


def main():
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("runs/F")
    csv_path = run_dir / "postProcessing/lagrangian/kinematicCloud/particles.csv"

    print("TEST 6 — live non-uniform air: particles vs a Python replay of the same steps")
    print(f"  run: {run_dir}")

    if not csv_path.is_file():
        raise SystemExit(f"cannot find the particle CSV: {csv_path}")

    # time folders written by the solver (0, 0.005, ... 1)
    times = sorted(
        (d.name for d in run_dir.iterdir()
         if d.is_dir() and d.name.replace(".", "", 1).isdigit()),
        key=float,
    )
    if len(times) < 2:
        raise SystemExit(f"need U in at least two time folders, found {times}")

    X, t_csv, d = load_particles(csv_path)
    n = X.shape[1]
    print(f"  U: {len(times)} time folders ({times[0]} .. {times[-1]} s), {NX*NY} cells")
    print(f"  particles: {n}, rows in CSV: {len(t_csv)}, t = {t_csv[0]:g} .. {t_csv[-1]:g} s")

    if len(times) != len(X):
        raise SystemExit(
            f"time folders ({len(times)}) and CSV rows per particle ({len(X)}) differ.\n"
            "The solver must write U and the particles at every step:\n"
            "  writeControl timeStep; writeInterval 1;"
        )
    if np.max(np.abs(d - D_P)) > 1e-12 * D_P:
        raise SystemExit(f"particle diameter {d} != case value {D_P}")

    U = [read_U(run_dir / s / "U") for s in times]
    changed = np.abs(U[-1] - U[0]).max()
    print(f"  air changed since t = 0: max |U(t_end) - U(0)| = {changed:.3e} m/s")
    if changed == 0.0:
        raise SystemExit(
            "the air never changed: this run cannot test the hand-over.\n"
            "Check that the lid moves and that U is written every step."
        )

    correct = replay(U, X)
    old_air = replay([U[max(s - 1, 0)] for s in range(len(U))], X)
    frozen = replay([U[0]] * len(U), X)

    e_correct = float(np.abs(X - correct).max())
    e_old = float(np.abs(X - old_air).max())
    e_frozen = float(np.abs(X - frozen).max())

    print()
    print(f"  {'Python uses':34s} {'max |x_NeoFOAM - x_Python|':>26s} "
          f"{'expected':>12s}  result")
    print(f"  {'-'*34:34s} {'-'*26:26s} {'-'*12:12s}  {'-'*6}")
    rows = [
        ("U(t_n)    this step (correct)", e_correct, f"< {TOL_X:g}", e_correct < TOL_X),
        ("U(t_n-1)  one step old", e_old, f"> {TOL_CONTROL:g}", e_old > TOL_CONTROL),
        ("U(0)      start air, frozen", e_frozen, f"> {TOL_CONTROL:g}",
         e_frozen > TOL_CONTROL),
    ]
    for name, err, expect, ok in rows:
        print(f"  {name:34s} {err:26.3e} {expect:>12s}  {'PASS' if ok else 'FAIL'}")

    moved = np.linalg.norm(X[-1] - X[0], axis=1)
    print()
    print(f"  distance moved by each particle [m]: "
          f"{np.array2string(moved, precision=4, floatmode='fixed')}")
    print(f"  particles still inside at t_end: {int(inside(X[-1]).sum())} of {n}")
    print(f"  max particle speed at t_end: "
          f"{np.linalg.norm(np.diff(X[-2:], axis=0)[0], axis=1).max() / DT:.3e} m/s")

    # The replay does not model the wall interaction (NeoFOAM rebounds a particle
    # that touches a wall). Report the margin so a run that grazes a wall is
    # recognisable instead of silently failing the comparison.
    clearance = float(min(
        X[:, :, 0].min(), LX - X[:, :, 0].max(),
        X[:, :, 1].min(), LY - X[:, :, 1].max(),
        X[:, :, 2].min(), LZ - X[:, :, 2].max(),
    ))
    print(f"  closest approach to a wall: {clearance:.4f} m (cell {DX:g} m)")
    if not inside(X.reshape(-1, 3)).all():
        print("  WARNING: a particle left the domain. NeoFOAM rebounds it")
        print("           (standardWallInteraction), the replay does not: the")
        print("           comparison is not meaningful for that particle.")

    print()
    if all(r[3] for r in rows):
        print("PASS: the particles follow the air of each step; a stale or frozen")
        print("      air is detected by the two controls.")
        return 0

    print("FAIL: see the table above.")
    print("  row 1 large  -> NeoFOAM did not use the air of this step")
    print("  row 2 small  -> the replay cannot detect one-step-old air")
    print("  row 3 small  -> the replay cannot detect frozen air")
    return 1


if __name__ == "__main__":
    sys.exit(main())
