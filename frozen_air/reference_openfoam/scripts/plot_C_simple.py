#!/usr/bin/env python3
"""Case C (sphereDrag + Euler): paths and velocity in x, y, z, NeoFOAM and OpenFOAM.

OpenFOAM moves an injected particle twice in its first step: its output at t is the NeoFOAM
state at t + dt, so OpenFOAM is plotted at t + dt (as compare.py, offset 1).

Writes figures/ref_C_simple.png. Run from the reference_openfoam folder:
    python3 scripts/plot_C_simple.py
"""
import numpy as np
import matplotlib.pyplot as plt

L = 1.0                                  # exit face x = 1 m
COLORS = ["tab:blue", "tab:orange", "tab:green"]


def read(code):
    D = np.genfromtxt(f"results/{code}/C.csv", delimiter=",", names=True)
    return {int(i): D[D["id"] == i] for i in np.unique(D["id"])}


N, O = read("neofoam"), read("openfoam")
dt = np.diff(np.unique(N[0]["t"]))[0]
for o in O.values():
    o["t"] += dt                         # the state OpenFOAM really holds (one step ahead)
fig = plt.figure(figsize=(16, 9))
gs = fig.add_gridspec(2, 6)
ax = [fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])]
av = [fig.add_subplot(gs[1, 2 * k:2 * k + 2]) for k in range(3)]

largest = 0.0
for pid, c in zip(sorted(N), COLORS):
    n, o = N[pid], O[pid]
    n = n[n["t"] <= o["t"][-1]]          # OpenFOAM removes a particle that leaves: plot both until then
    lbl = f"P{pid}, d = {n['d'][0] * 1e6:.0f} µm"
    for a, (i, j) in zip(ax, [("x", "z"), ("x", "y")]):
        a.plot(n[i], n[j], "-", color=c, lw=2, label=lbl)
        a.plot(o[i], o[j], "o", color=c, ms=4, mfc="none", markevery=20)
    for a, k in zip(av, ["ux", "uy", "uz"]):
        a.plot(n["t"], n[k], "-", color=c, lw=2, label=lbl)
        a.plot(o["t"], o[k], "o", color=c, ms=4, mfc="none", markevery=20)
    # speed difference at the same output times
    common, i_n, i_o = np.intersect1d(np.round(n["t"], 9), np.round(o["t"], 9), return_indices=True)
    for k in ["ux", "uy", "uz"]:
        largest = max(largest, np.abs(n[k][i_n] - o[k][i_o]).max())

for a in ax:
    a.axvline(L, color="k", lw=1)
    a.text(L, a.get_ylim()[0], " exit face x = 1 m", rotation=90, va="bottom", fontsize=8)
ax[0].set(xlabel="x [m]", ylabel="z [m]", title="Side view (x–z): wind →, gravity ↓")
ax[1].set(xlabel="x [m]", ylabel="y [m]", title="Top view (x–y)")
for a, n in zip(av, "xyz"):
    a.set(xlabel="t [s]", ylabel=f"u_{n} [m/s]", title=f"Velocity in {n}")
for a in (ax[0], av[0]):
    a.plot([], [], "k-", lw=2, label="NeoFOAM (lines)")
    a.plot([], [], "ko", ms=4, mfc="none", label="OpenFOAM (circles)")
    a.legend(fontsize=8)
fig.suptitle(f"Case C: three particles in the wind, 10×10×10 cells (sphereDrag + Euler) — "
             f"NeoFOAM and OpenFOAM velocities differ by at most {largest:.1e} m/s", fontsize=12)
fig.tight_layout()
fig.savefig("figures/ref_C_simple.png", dpi=150)
print(f"largest velocity difference: {largest:.2e} m/s")
print("written figures/ref_C_simple.png")
