#!/usr/bin/env python3
"""Compare NeoFOAM with the OpenFOAM v2406 reference for test 6 (sphereDrag + Euler, tolerance 1e-12).

1. Air: NeoFOAM (runs/neofoam/F) vs OpenFOAM icoFoam (runs/openfoam/F_fluid), every step,
   max |U_NeoFOAM - U_OpenFOAM| over all cells.
2. Particles in the same air: NeoFOAM vs OpenFOAM icoUncoupledKinematicParcelFoam driven step by
   step with the NeoFOAM air (results/openfoam/F_chain.csv). OpenFOAM at t_k+1 is the NeoFOAM
   state at t_k (the first OpenFOAM step runs in the start air U(0) = 0).
   OpenFOAM splits a step where a particle crosses a face and uses the air of each cell for its
   part of the step; NeoFOAM takes the air of the start cell for the whole step. So the two must
   agree to round-off until a particle first crosses a face, and differ after.
3. Information only: OpenFOAM's own one-way solver kinematicParcelFoam (results/openfoam/F.csv):
   its air is pimpleFoam's (not icoFoam's) and it moves the particles before solving the air.

Run from the reference_openfoam folder:   python3 scripts/compare.py | tee results/summary.txt
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from read_openfoam import read_U  # noqa: E402

DT = 5e-3
DX = 0.05                                 # cell size, 20 x 20 cells
TOL_FLUID = 1e-9                          # m/s, same equations at tolerance 1e-12
TOL_SAME = 1e-12                          # m, particles before the first face crossing


def times(run):
    return sorted((d for d in Path(run).iterdir() if d.is_dir() and (d / "U").is_file()
                   and d.name.replace(".", "", 1).replace("e-", "", 1).isdigit()),
                  key=lambda d: float(d.name))


def air(run):
    """{time step index: cell values of U [400, 3]}"""
    out = {}
    for d in times(run):
        u = read_U(d / "U")
        out[round(float(d.name) / DT)] = np.zeros((400, 3)) if u is None else np.array(u)
    return out


def particles(path):
    D = np.genfromtxt(path, delimiter=",", names=True)
    return {int(i): D[D["id"] == i] for i in np.unique(D["id"])}


def align(n, o, shift):
    """indices of the common states: NeoFOAM step k  <->  OpenFOAM step k + shift"""
    kn = np.round(n["t"] / DT).astype(int)
    ko = np.round(o["t"] / DT).astype(int) - shift
    return np.intersect1d(kn, ko, return_indices=True)


def first_crossing(n):
    """time step index of the first face crossing (cell of the particle changes)"""
    cell = np.floor(n["x"] / DX) + 20 * np.floor(n["y"] / DX)
    k = np.flatnonzero(np.diff(cell) != 0)
    return int(round(n["t"][k[0] + 1] / DT)) if len(k) else None


def main():
    ok = True
    print("TEST 6 reference: NeoFOAM vs OpenFOAM v2406, lid-driven cavity, sphereDrag + Euler, tolerance 1e-12\n")

    # 1. air
    N, O = air("runs/neofoam/F"), air("runs/openfoam/F_fluid")
    common = sorted(set(N) & set(O))
    dU = np.array([np.abs(N[k] - O[k]).max() for k in common])
    k_end = common[-1]
    verdict = "PASS" if dU.max() < TOL_FLUID else "FAIL"
    ok &= verdict == "PASS"
    print("1. AIR  NeoFOAM vs OpenFOAM icoFoam, every step")
    print(f"   steps compared: {len(common)}, max |U| at t_end: {np.linalg.norm(O[k_end], axis=1).max():.3f} m/s")
    print(f"   max |U_NeoFOAM - U_OpenFOAM| over all cells and steps: {dU.max():.2e} m/s   "
          f"(< {TOL_FLUID:g})  {verdict}\n")

    # 2. particles in the same air
    Np, Op = particles("results/neofoam/F.csv"), particles("results/openfoam/F_chain.csv")
    print("2. PARTICLES in the same air: NeoFOAM vs OpenFOAM icoUncoupledKinematicParcelFoam (NeoFOAM U every step)")
    print("   P   first face crossing   max|dx| before it [m]   max|dx| after [m]   max|du| after [m/s]   verdict")
    for p in sorted(Np):
        n, o = Np[p], Op[p]
        c, i_n, i_o = align(n, o, 1)
        dx = np.max([np.abs(n[k][i_n] - o[k][i_o]) for k in "xyz"], axis=0)
        du = np.max([np.abs(n[k][i_n] - o[k][i_o]) for k in ("ux", "uy", "uz")], axis=0)
        kc = first_crossing(n)
        before = c < kc if kc is not None else np.ones_like(c, bool)
        v = "PASS" if dx[before].max() < TOL_SAME else "FAIL"
        ok &= v == "PASS"
        after_x = f"{dx[~before].max():.2e}" if (~before).any() else "-"
        after_u = f"{du[~before].max():.2e}" if (~before).any() else "-"
        tc = f"t = {kc * DT:.3f} s" if kc is not None else "none"
        print(f"   P{p}  {tc:>20s}   {dx[before].max():>21.2e}   {after_x:>17s}   {after_u:>19s}   {v}")
    print(f"   (pass: identical to {TOL_SAME:g} m until the first face crossing; after it OpenFOAM splits the step"
          " at the face, so a difference is expected)\n")

    # 3. information: OpenFOAM's own one-way solver
    Kp = particles("results/openfoam/F.csv")
    K = air("runs/openfoam/F")
    dK = max(np.abs(K[k] - O[k]).max() for k in common if k in K)
    worst = 0.0
    for p in sorted(Np):
        c, i_n, i_o = align(Np[p], Kp[p], 1)
        worst = max(worst, max(np.abs(Np[p][k][i_n] - Kp[p][k][i_o]).max() for k in "xyz"))
    print("3. INFORMATION: OpenFOAM kinematicParcelFoam (its own one-way solver)")
    print(f"   its air vs icoFoam: max |dU| = {dK:.2e} m/s (pimpleFoam momentum equation, not icoFoam's)")
    print(f"   its particles vs NeoFOAM: max |dx| = {worst:.2e} m (different air + particles moved before the air is solved)\n")

    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
