#!/usr/bin/env python3
"""Case F, air: NeoFOAM vs OpenFOAM icoFoam (same cavity, tolerance 1e-12).

Left, middle: velocity along the centre lines of the cavity at three times (the classic cavity
plot): u_x along the vertical line x = 0.5 m, u_y along the horizontal line y = 0.5 m. The centre
lines lie on cell faces (20 cells), so the two neighbouring cell columns/rows are averaged.
Right: largest difference over all cells at every step.

Writes figures/ref_F_fluid.png. Run from the reference_openfoam folder:
    python3 scripts/plot_F_fluid.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from compare import DT, air  # noqa: E402

NX = 20
TIMES = [0.1, 0.5, 3.0]
COLORS = ["tab:blue", "tab:orange", "tab:green"]

N, O = air("runs/neofoam/F"), air("runs/openfoam/F_fluid")
common = sorted(set(N) & set(O))
centres = (np.arange(NX) + 0.5) / NX


def centrelines(U):
    G = U.reshape(NX, NX, 3)                       # [j (y), i (x), component], cell = i + 20 j
    ux = 0.5 * (G[:, NX // 2 - 1, 0] + G[:, NX // 2, 0])     # along y at x = 0.5
    uy = 0.5 * (G[NX // 2 - 1, :, 1] + G[NX // 2, :, 1])     # along x at y = 0.5
    return ux, uy


fig, ax = plt.subplots(1, 3, figsize=(16, 5))
for t, c in zip(TIMES, COLORS):
    k = round(t / DT)
    (nx_, ny_), (ox_, oy_) = centrelines(N[k]), centrelines(O[k])
    ax[0].plot(nx_, centres, "-", color=c, lw=2, label=f"t = {t:g} s")
    ax[0].plot(ox_, centres, "o", color=c, ms=5, mfc="none")
    ax[1].plot(centres, ny_, "-", color=c, lw=2, label=f"t = {t:g} s")
    ax[1].plot(centres, oy_, "o", color=c, ms=5, mfc="none")

dU = np.array([np.abs(N[k] - O[k]).max() for k in common])
t = np.array(common) * DT
ax[2].semilogy(t[1:], dU[1:], "k-", lw=1)
ax[2].axhline(1e-9, color="k", ls="--", lw=1, label="pass limit 1e-9 m/s")
ax[2].set_ylim(1e-14, 1e-6)

ax[0].set(xlabel="u_x [m/s]", ylabel="y [m]", title="u_x along the vertical centre line x = 0.5 m")
ax[1].set(xlabel="x [m]", ylabel="u_y [m/s]", title="u_y along the horizontal centre line y = 0.5 m")
ax[2].set(xlabel="t [s]", ylabel="max |U_NeoFOAM − U_OpenFOAM| [m/s]",
          title="Largest difference over all 400 cells, every step")
for a in ax[:2]:
    a.axhline(0, color="0.7", lw=0.8) if a is ax[1] else a.axvline(0, color="0.7", lw=0.8)
    a.plot([], [], "k-", lw=2, label="NeoFOAM (lines)")
    a.plot([], [], "ko", ms=5, mfc="none", label="OpenFOAM icoFoam (circles)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3, which="both")
fig.suptitle(f"Case F, air: lid-driven cavity, Re = 100, air starts at rest — NeoFOAM and OpenFOAM icoFoam "
             f"differ by at most {dU.max():.1e} m/s", fontsize=12)
fig.tight_layout()
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/ref_F_fluid.png", dpi=150)
print(f"largest air difference: {dU.max():.2e} m/s")
print("written figures/ref_F_fluid.png")
