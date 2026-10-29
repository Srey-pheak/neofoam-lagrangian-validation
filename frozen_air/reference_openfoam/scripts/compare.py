#!/usr/bin/env python3
"""Compare NeoFOAM with the OpenFOAM v2406 reference (same cases: sphereDrag + Euler).

Reads results/openfoam/<name>.csv and results/neofoam/<name>.csv (t,id,x,y,z,ux,uy,uz,d),
matches the particles by diameter and compares them at the same physical states.

OpenFOAM v2406 moves a particle injected at SOI 0 twice in its first time step, so its
output at t_k is expected to equal the NeoFOAM state one step later. Both alignments
(offset 0 and 1 step) are compared and the better one is reported.

Also, per code:
  A  error vs a fine RK4 reference of the same equation (Putnam drag, no gravity) and
     the order of accuracy over the dt sweep
  B  terminal velocity at t_end vs the force balance
  C  exit time of the particle that leaves through x = 1 m

Run from the reference_openfoam folder:   python3 scripts/compare.py
Writes figures/ref_*.png; pipe through tee for results/summary.txt.
"""
import glob
import math
import os
import re

import numpy as np
import matplotlib.pyplot as plt

rho_p, rho_f, nu, g = 2500.0, 1.2, 1.5e-5, 9.81
mu = rho_f * nu
keys = ["x", "y", "z", "ux", "uy", "uz"]
os.makedirs("figures", exist_ok=True)


def CdRe(Re):
    return 0.424 * Re if Re > 1000 else 24 * (1 + Re ** (2 / 3) / 6)


def read(path):
    D = np.genfromtxt(path, delimiter=",", names=True)
    return {dd: D[D["d"] == dd] for dd in np.unique(D["d"])}


def compare(n, o, dt):
    """max relative velocity difference for offset 0 and 1 step -> (best offset, diff, idx n, idx o)"""
    best = None
    for off in (0, 1):
        tn = np.round(n["t"] / dt).astype(int)
        to = np.round(o["t"] / dt).astype(int) + off
        common, i_n, i_o = np.intersect1d(tn, to, return_indices=True)
        if len(common) == 0:
            continue
        uref = max(np.abs(o[k][i_o]).max() for k in keys[3:])
        du = max(np.abs(n[k][i_n] - o[k][i_o]).max() for k in keys[3:]) / uref
        dx = max(np.abs(n[k][i_n] - o[k][i_o]).max() for k in keys[:3])
        if best is None or du < best[1]:
            best = (off, du, dx, i_n, i_o)
    return best


def rk4_relaxation(d, u0, t_end, dt=2e-6):
    """reference u(t) for test 1: Putnam drag, air at rest, no gravity"""
    m = rho_p * math.pi * d**3 / 6

    def a(u):
        Re = rho_f * d * abs(u) / mu
        return -math.pi / 8 * mu * d * CdRe(Re) * u / m

    n = int(round(t_end / dt))
    ts, us, u = [0.0], [u0], u0
    for k in range(n):
        k1 = a(u)
        k2 = a(u + dt / 2 * k1)
        k3 = a(u + dt / 2 * k2)
        k4 = a(u + dt * k3)
        u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ts.append((k + 1) * dt)
        us.append(u)
    return np.array(ts), np.array(us)


def terminal_velocity(d):
    vs = rho_p * d**2 / (18 * mu) * g * (1 - rho_f / rho_p)
    v = vs
    for _ in range(1000):
        v = vs * 24 / CdRe(rho_f * v * d / mu)
    return v


names = sorted(os.path.basename(f)[:-4] for f in glob.glob("results/openfoam/*.csv"))
names = [n for n in names if os.path.exists(f"results/neofoam/{n}.csv")]
if not names:
    raise SystemExit("no matching results: run scripts/run_reference.sh first")

# ─── 1. NeoFOAM vs OpenFOAM, every case ───
print("1. NeoFOAM vs OpenFOAM (same case, sphereDrag + Euler)")
print(f"{'case':12s} {'d [um]':>7s} {'offset':>7s} {'max|du|/max|u|':>15s} {'max|dx| [m]':>12s}")
parity = {}
for name in names:
    N, O = read(f"results/neofoam/{name}.csv"), read(f"results/openfoam/{name}.csv")
    dt = np.diff(np.unique(np.concatenate([v["t"] for v in N.values()])))[0]
    for dd in sorted(N):
        if dd not in O:
            print(f"{name:12s} {dd * 1e6:7.0f}  missing in the OpenFOAM output")
            continue
        off, du, dx, i_n, i_o = compare(N[dd], O[dd], dt)
        parity[(name, dd)] = (off, du, i_n, i_o, dt)
        print(f"{name:12s} {dd * 1e6:7.0f} {off:7d} {du:15.3e} {dx:12.3e}")
print("""notes: offset 1 = OpenFOAM one step ahead (double first step of an injected particle).
       B is written every 40 steps, so the 1-step offset cannot be applied there: its difference
       is the transient seen one step apart; compare its terminal velocities (section 3).
       C crosses faces: OpenFOAM splits a step at every face (maxCo), NeoFOAM does not.""")

# ─── 2. A: error vs fine reference and order, both codes ───
A = sorted((float(re.search(r"dt([0-9.eE+-]+)", n).group(1)), n) for n in names if n.startswith("A_"))
if A:
    print("\n2. TEST 1 CASE with sphereDrag: error vs fine RK4 reference, both codes")
    d0, u0 = 5e-5, 0.1
    tr, ur = rk4_relaxation(d0, u0, 0.2)
    print(f"{'dt [s]':>10s} {'e NeoFOAM':>11s} {'p':>5s} {'e OpenFOAM':>11s} {'p':>5s}")
    rows, prev = [], None
    for dt, name in sorted(A, reverse=True):
        es = []
        for code in ("neofoam", "openfoam"):
            D = read(f"results/{code}/{name}.csv")[d0]
            # OpenFOAM: compare at the state it really holds (one step ahead), see section 1
            shift = parity[(name, d0)][0] * dt if code == "openfoam" else 0.0
            ue = np.interp(D["t"] + shift, tr, ur)
            es.append(np.abs(D["ux"] - ue).max() / u0)
        p = ["", ""] if prev is None else [f"{math.log(prev[1][k] / es[k]) / math.log(prev[0] / dt):.2f}" for k in (0, 1)]
        print(f"{dt:10.3e} {es[0]:11.3e} {p[0]:>5s} {es[1]:11.3e} {p[1]:>5s}")
        rows.append((dt, *es))
        prev = (dt, es)
    r = np.array(rows)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(r[:, 0], r[:, 1], "o-", ms=7, label="NeoFOAM Euler")
    ax.loglog(r[:, 0], r[:, 2], "s--", ms=9, mfc="none", label="OpenFOAM Euler")
    ax.loglog(r[:, 0], r[-1, 1] * (r[:, 0] / r[-1, 0]), "k:", lw=1, label="slope 1")
    ax.set(xlabel="dt [s]", ylabel="error = max $|u - u_{reference}|\\ /\\ u_0$", title="Test 1 case with sphereDrag: error of both codes")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/ref_A_convergence.png", dpi=150)

# ─── 3. B: terminal velocity, both codes ───
if "B" in names:
    print("\n3. TEST 2 CASE: terminal velocity at t_end vs force balance, both codes")
    print(f"{'d [um]':>7s} {'v_t ref':>12s} {'NeoFOAM':>12s} {'err':>9s} {'OpenFOAM':>12s} {'err':>9s}")
    N, O = read("results/neofoam/B.csv"), read("results/openfoam/B.csv")
    for dd in sorted(N):
        vt = terminal_velocity(dd)
        vn, vo = -N[dd]["uz"][-1], -O[dd]["uz"][-1]
        print(f"{dd * 1e6:7.0f} {vt:12.6e} {vn:12.6e} {abs(vn - vt) / vt:9.1e} {vo:12.6e} {abs(vo - vt) / vt:9.1e}")

# ─── 4. C: exit time of the particle that leaves ───
if "C" in names:
    print("\n4. TEST 3 CASE: particles in the wind (crossing cells), exit through x = 1 m")
    N, O = read("results/neofoam/C.csv"), read("results/openfoam/C.csv")
    for dd in sorted(N):
        n = N[dd]
        kn = np.argmax(n["x"] >= 1.0) if (n["x"] >= 1.0).any() else None
        o = O.get(dd)
        lo = o["t"][-1] if o is not None else float("nan")
        if kn is not None:
            print(f"  d = {dd * 1e6:.0f} um: NeoFOAM first x >= 1 m at t = {n['t'][kn]:.4f} s; "
                  f"OpenFOAM last output with the particle t = {lo:.4f} s (removed when it leaves)")
        else:
            print(f"  d = {dd * 1e6:.0f} um: stays inside (NeoFOAM end t = {n['t'][-1]:.3f} s, OpenFOAM {lo:.3f} s)")

print("\nfigures: figures/ref_A_convergence.png")
