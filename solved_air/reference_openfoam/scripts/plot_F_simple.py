#!/usr/bin/env python3
"""Case F, particles in the same air: NeoFOAM vs OpenFOAM icoUncoupledKinematicParcelFoam driven
step by step with the NeoFOAM air (sphereDrag + Euler).

OpenFOAM at t_k+1 is the NeoFOAM state at t_k (its first step runs in the start air), so OpenFOAM
is plotted at t - dt. Both codes agree to round-off until a particle first crosses a cell face;
there OpenFOAM splits the step and uses the air of each cell for its part, NeoFOAM uses the air of
the start cell for the whole step.

Writes figures/ref_F_simple.png. Run from the reference_openfoam folder:
    python3 scripts/plot_F_simple.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from compare import DT, align, first_crossing, particles  # noqa: E402

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

N, O = particles("results/neofoam/F.csv"), particles("results/openfoam/F_chain.csv")
fig = plt.figure(figsize=(16, 9))
gs = fig.add_gridspec(2, 2)
ap, ae = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
gv = gs[1, :].subgridspec(1, 2)
av = [fig.add_subplot(gv[0, 0]), fig.add_subplot(gv[0, 1])]

largest = 0.0
for pid, c in zip(sorted(N), COLORS):
    n, o = N[pid], O[pid].copy()
    o["t"] -= DT                                   # the NeoFOAM state OpenFOAM holds
    lbl = f"P{pid}"
    ap.plot(n["x"], n["y"], "-", color=c, lw=2, label=lbl)
    ap.plot(o["x"], o["y"], "o", color=c, ms=4, mfc="none", markevery=20)
    ap.plot(n["x"][0], n["y"][0], "o", color=c, ms=7)
    for a, k in zip(av, ["ux", "uy"]):
        a.plot(n["t"], n[k], "-", color=c, lw=1.5, label=lbl)
        a.plot(o["t"], o[k], "o", color=c, ms=3, mfc="none", markevery=15)
    common, i_n, i_o = align(N[pid], O[pid], 1)
    dx = np.max([np.abs(N[pid][k][i_n] - O[pid][k][i_o]) for k in "xyz"], axis=0)
    largest = max(largest, dx.max())
    ae.semilogy(common * DT, np.maximum(dx, 1e-18), "-", color=c, lw=1.2, label=lbl)
    kc = first_crossing(N[pid])
    ae.axvline(kc * DT, color=c, ls=":", lw=1)

ap.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, lw=2))
ap.plot([0, 1], [1, 1], color="tab:red", lw=4)
ap.text(0.5, 1.02, "lid U = (1, 0, 0) m/s", color="tab:red", ha="center", fontsize=9)
for g in np.arange(0.05, 1.0, 0.05):
    ap.axvline(g, color="0.9", lw=0.5, zorder=0)
    ap.axhline(g, color="0.9", lw=0.5, zorder=0)
ap.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.08), aspect="equal", xlabel="x [m]", ylabel="y [m]",
       title="Paths, t = 0 … 3 s (● start, grey: cell faces)")
ae.axhline(1e-12, color="k", ls="--", lw=1, label="pass limit 1e-12 m (before the first face crossing)")
ae.plot([], [], "k:", lw=1, label="first face crossing of each particle")
ae.set_ylim(1e-18, 1)
ae.set(xlabel="t [s]", ylabel="|x_NeoFOAM − x_OpenFOAM| [m]",
       title="Position difference: round-off until the first face crossing")
ae.legend(fontsize=8, loc="lower right")
ae.grid(alpha=0.3, which="both")
for a, k in zip(av, "xy"):
    a.set(xlabel="t [s]", ylabel=f"u_{k} [m/s]", title=f"Particle velocity in {k}")
    a.grid(alpha=0.3)
for a in (ap, av[0]):
    a.plot([], [], "k-", lw=2, label="NeoFOAM (lines)")
    a.plot([], [], "ko", ms=4, mfc="none", label="OpenFOAM (circles)")
    a.legend(fontsize=8, ncol=2)
fig.suptitle("Case F, particles: 5 particles (d = 1 mm) in the same solved cavity air (sphereDrag + Euler) — "
             f"identical until a particle crosses a face, at most {largest:.1e} m apart after", fontsize=12)
fig.tight_layout()
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/ref_F_simple.png", dpi=150)
print(f"largest position difference: {largest:.2e} m")
print("written figures/ref_F_simple.png")
