#!/usr/bin/env python3
"""Case B (sphereDrag + Euler): falling speed over time and final speed, NeoFOAM and OpenFOAM.

Writes figures/ref_B_simple.png. Run from the reference_openfoam folder:
    python3 scripts/plot_B_simple.py
"""
import numpy as np
import matplotlib.pyplot as plt

rho_p, rho_f, nu, g = 2500.0, 1.2, 1.5e-5, 9.81
mu = rho_f * nu
g_eff = g * (1 - rho_f / rho_p)
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]


def CdRe(Re):
    return 0.424 * Re if Re > 1000 else 24 * (1 + Re ** (2 / 3) / 6)


def final_speed(d):
    """theory: gravity = drag (Putnam), solved by fixed-point iteration"""
    vs = rho_p * d**2 / (18 * mu) * g_eff
    v = vs
    for _ in range(1000):
        v = vs * 24 / CdRe(rho_f * v * d / mu)
    return v


def read(code):
    D = np.genfromtxt(f"results/{code}/B.csv", delimiter=",", names=True)
    return {dd: D[D["d"] == dd] for dd in np.unique(D["d"])}


N, O = read("neofoam"), read("openfoam")
sizes = sorted(N)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

# 1. falling speed over time, both codes
for dd, c in zip(sizes, COLORS):
    ax[0].plot(N[dd]["t"], -N[dd]["uz"], "-", color=c, lw=2, label=f"d = {dd * 1e6:.0f} µm")
    ax[0].plot(O[dd]["t"], -O[dd]["uz"], "o", color=c, ms=4, mfc="none", markevery=25)
ax[0].plot([], [], "k-", lw=2, label="NeoFOAM (lines)")
ax[0].plot([], [], "ko", ms=4, mfc="none", label="OpenFOAM (circles)")
ax[0].set(xlabel="t [s]", ylabel="falling speed [m/s]", xlim=(0, 2.5),
          title="Falling speed over time: both codes on top of each other")
ax[0].legend(fontsize=8, loc="center right")

# 2. final speed of each size: theory, NeoFOAM, OpenFOAM (bars side by side, numbers on top)
x = np.arange(len(sizes))
w = 0.27
theory = [final_speed(dd) for dd in sizes]
neo = [-N[dd]["uz"][-1] for dd in sizes]
of = [-O[dd]["uz"][-1] for dd in sizes]
for k, (vals, lbl, col) in enumerate([(theory, "theory (gravity = drag)", "lightgray"),
                                      (neo, "NeoFOAM", "tab:blue"),
                                      (of, "OpenFOAM", "tab:orange")]):
    bars = ax[1].bar(x + (k - 1) * w, vals, w, label=lbl, color=col, edgecolor="k", lw=0.5)
    if k == 0:
        for b, v in zip(bars, vals):
            ax[1].text(b.get_x() + 1.5 * w, v * 1.06, f"{v:.3g} m/s", ha="center", fontsize=8)
ax[1].set_yscale("log")
ax[1].set_xticks(x, [f"{dd * 1e6:.0f} µm" for dd in sizes])
ax[1].set(xlabel="particle diameter", ylabel="final falling speed [m/s]  (log scale)",
          title="Final speed: theory, NeoFOAM and OpenFOAM are the same", ylim=(4e-3, 8))
ax[1].legend(fontsize=8, loc="upper left")

largest = max(abs(a - b) / b for a, b in zip(neo, of))
fig.suptitle(f"Case B: five particles falling from rest (sphereDrag + Euler) — "
             f"NeoFOAM and OpenFOAM final speeds differ by at most {largest:.0e}", fontsize=12)
fig.tight_layout()
fig.savefig("figures/ref_B_simple.png", dpi=150)
for dd, t, n, o in zip(sizes, theory, neo, of):
    print(f"d = {dd * 1e6:3.0f} um   theory {t:.6f}   NeoFOAM {n:.6f}   OpenFOAM {o:.6f}  m/s")
print("written figures/ref_B_simple.png")
