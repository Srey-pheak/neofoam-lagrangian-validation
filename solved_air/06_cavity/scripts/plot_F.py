#!/usr/bin/env python3
"""Test 6 — figures: particle paths in the solved cavity air, and the controls.

Left: the paths of the 5 particles over the cavity at t_end, with the solved air
(speed + direction per cell) as the background. This is the picture of "the
particles follow the air of the cavity".

Right: |x_NeoFOAM - x_Python| against time for the three replays of
scripts/replay_F.py. The correct air stays at round-off; the two wrong-air
controls (one step old, frozen at t = 0) grow. That is the visual proof that the
test can detect a stale U.

Writes figures/F_particles.png.

Usage:
    python3 scripts/plot_F.py [runs/06_cavity]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_F as R  # noqa: E402

COL = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "runs/06_cavity"
    csv_path = os.path.join(run_dir, "postProcessing/lagrangian/kinematicCloud/particles.csv")
    if not os.path.isdir(run_dir) or not os.path.isfile(csv_path):
        sys.exit(f"no run in {run_dir}: run scripts/run_F.sh first")

    times = sorted(
        (d for d in os.listdir(run_dir)
         if os.path.isdir(os.path.join(run_dir, d)) and d.replace(".", "", 1).isdigit()),
        key=float,
    )
    if len(times) < 2:
        sys.exit(f"need U in at least two time folders in {run_dir}")

    X, t, d = R.load_particles(csv_path)
    U = [R.read_U(os.path.join(run_dir, s, "U")) for s in times]
    U_end = U[-1]

    # ─── figure ───
    fig = plt.figure(figsize=(13, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.42)
    ax = fig.add_subplot(gs[0, 0])

    # background: speed of the solved air on the cell centres
    i_c, j_c = np.meshgrid(np.arange(R.NX), np.arange(R.NY))
    speed = np.linalg.norm(U_end, axis=1).reshape(R.NY, R.NX)
    xc = (i_c + 0.5) * R.DX
    yc = (j_c + 0.5) * R.DX
    pm = ax.pcolormesh(
        np.arange(R.NX + 1) * R.DX, np.arange(R.NY + 1) * R.DX, speed,
        cmap="Blues", shading="flat", alpha=0.85,
    )
    cb = fig.colorbar(pm, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.set_title("|U| [m/s]", fontsize=8.5, pad=6)

    # air direction per cell
    ax.quiver(
        xc, yc, U_end[:, 0].reshape(R.NY, R.NX), U_end[:, 1].reshape(R.NY, R.NX),
        color="0.25", alpha=0.75, width=0.003, scale=18,
    )

    # particle paths
    for p in range(X.shape[1]):
        c = COL[p % len(COL)]
        ax.plot(X[:, p, 0], X[:, p, 1], color=c, lw=1.7, zorder=4,
                label=f"P{p}  d = {d[p]*1e3:.1f} mm")
        ax.plot(X[0, p, 0], X[0, p, 1], "o", color=c, ms=8, mec="white", mew=1.0, zorder=5)
        ax.plot(X[-1, p, 0], X[-1, p, 1], "^", color=c, ms=9, mec="white", mew=1.0, zorder=5)

    # cavity outline + lid
    ax.add_patch(Rectangle((0, 0), R.LX, R.LY, fill=False, ec="k", lw=2.0, zorder=6))
    ax.plot([0, 1], [1, 1], color="tab:red", lw=4, solid_capstyle="butt", zorder=7)
    ax.annotate("", xy=(0.62, 1.045), xytext=(0.38, 1.045),
                arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2.0))
    ax.text(0.5, 1.075, "lid U = (1, 0, 0) m/s", color="tab:red", fontsize=9.5,
            ha="center", va="bottom")

    ax.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.12), aspect="equal",
           xlabel="x [m]", ylabel="y [m]")
    ax.set_title(f"Test 6 — particle paths, {os.path.basename(run_dir)}\n"
                 f"t = {t[0]:g} … {t[-1]:g} s, ○ start, ▲ end", fontsize=11)
    ax.legend(loc="lower left", fontsize=8.5, framealpha=0.9)

    # ─── right: the three replays ───
    ax2 = fig.add_subplot(gs[0, 1])
    replays = {
        "U(t_n) — this step (correct)": (R.replay(U, X), "tab:green", "-"),
        "U(t_n-1) — one step old": (R.replay([U[max(s - 1, 0)] for s in range(len(U))], X),
                                     "tab:orange", "--"),
        "U(0) — start air (frozen)": (R.replay([U[0]] * len(U), X), "tab:red", ":"),
    }
    for label, (ref, c, ls) in replays.items():
        err = np.abs(X - ref).max(axis=(1, 2))     # worst particle, per step
        ax2.semilogy(t, np.maximum(err, 1e-18), color=c, ls=ls, lw=1.8, label=label)

    ax2.axhline(1e-9, color="k", lw=1.0, ls="--", label="pass limit 1e-9 m")
    ax2.axhline(1e-6, color="0.5", lw=1.0, ls=":", label="control limit 1e-6 m")
    ax2.set(xlabel="t [s]", ylabel="max |Δx| [m]",
            ylim=(1e-18, 2.0), xlim=(t[0], t[-1]))
    ax2.set_title("The three replays of the same steps\n"
                  "correct air: round-off; wrong air: detected", fontsize=11)
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(loc="lower right", fontsize=8.5)

    out = os.path.join("figures", "F_particles.png")
    os.makedirs("figures", exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  wrote {out}  ({len(times)} time folders, {X.shape[1]} particles)")

    print(f"  air speed at t_end: max {speed.max():.3f} m/s, "
          f"min {speed.min():.3f} m/s, mean {speed.mean():.3f} m/s")
    for label, (ref, _, _) in replays.items():
        print(f"  {label:32s} max |Δx| = {np.abs(X - ref).max():.3e} m")


if __name__ == "__main__":
    main()
