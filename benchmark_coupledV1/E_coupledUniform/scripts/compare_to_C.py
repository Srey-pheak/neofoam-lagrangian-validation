#!/usr/bin/env python3
"""Test E — coupled uniform flow vs Test C (frozen uniform air).

Test C moves the particles in air frozen at (2 0 0) m/s.  Test E solves the
same air every step (PISO + particles, one-way) and moves the particles with
the solved air.  A uniform flow is an exact solution of the incompressible
equations, so the particles must follow Test C's paths up to solver round-off.

Pass: same rows, same t and id columns, and
      max |x_E - x_C| < 1e-9 m, max |u_E - u_C| < 1e-9 m/s, d identical.

Usage:
    python3 scripts/compare_to_C.py [E.csv] [C.csv]

Defaults:
    results/E_uniform.csv
    ../../benchmark_uncoupleV1/C_3dBoundaries/results/C_RK2.csv
"""

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parents[1]   # .../E_coupledUniform
BENCH = HERE.parents[1]                       # .../Benchmarks

E_CANDIDATES = [
    HERE / "results/E_uniform.csv",
    HERE / "runs/E_uniform/postProcessing/lagrangian/kinematicCloud/particles.csv",
]
C_CANDIDATES = [
    BENCH / "benchmark_uncoupleV1/C_3dBoundaries/results/C_RK2.csv",
    BENCH / "C_3dBoundaries/results/C_RK2.csv",  # old layout
]

TOL_X = 1e-9  # m
TOL_U = 1e-9  # m/s


def first_existing(paths, what):
    for p in paths:
        if p.is_file():
            return p
    raise SystemExit(
        f"cannot find the {what} CSV, tried:\n  " + "\n  ".join(map(str, paths))
    )


def load(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    e_path = Path(sys.argv[1]) if len(sys.argv) > 1 else first_existing(E_CANDIDATES, "test E")
    c_path = Path(sys.argv[2]) if len(sys.argv) > 2 else first_existing(C_CANDIDATES, "test C")

    print("TEST E — one-way coupling, uniform flow (PISO + particles) vs Test C")
    print(f"  E: {e_path}")
    print(f"  C: {c_path}")

    E = load(e_path)
    C = load(c_path)
    print(f"  rows: E = {len(E)}, C = {len(C)}")
    print(f"  t:    E = {E['t'][0]:g} .. {E['t'][-1]:g} s, C = {C['t'][0]:g} .. {C['t'][-1]:g} s")

    ok = len(E) == len(C)
    if not ok:
        print("  FAIL: different number of rows")
    elif E.dtype.names != C.dtype.names:
        print(f"  FAIL: different columns: {E.dtype.names} vs {C.dtype.names}")
        ok = False
    else:
        print(f"  particles: {sorted({int(i) for i in E['id']})}")
        if not np.array_equal(E["t"], C["t"]):
            print("  FAIL: time column differs")
            ok = False
        if not np.array_equal(E["id"], C["id"]):
            print("  FAIL: particle id column differs")
            ok = False

        print()
        print(f"  {'column':>6} {'max |E-C|':>12} {'limit':>10}  result")
        print(f"  {'-'*6:>6} {'-'*12:>12} {'-'*10:>10}  {'-'*6}")
        for col, tol in (
            ("x", TOL_X), ("y", TOL_X), ("z", TOL_X),
            ("ux", TOL_U), ("uy", TOL_U), ("uz", TOL_U),
            ("d", 0.0),
        ):
            worst = float(np.abs(E[col] - C[col]).max()) if len(E) else 0.0
            passed = worst <= tol
            ok = ok and passed
            print(f"  {col:>6} {worst:12.3e} {tol:10.1e}  {'PASS' if passed else 'FAIL'}")

    print()
    if ok:
        print("PASS: coupled uniform flow reproduces Test C within round-off.")
        print("NOTE: a steady uniform flow cannot prove that U is re-read every step;")
        print("      a non-uniform or time-changing flow is the next test.")
        return 0

    print("FAIL: coupled uniform flow differs from Test C.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
