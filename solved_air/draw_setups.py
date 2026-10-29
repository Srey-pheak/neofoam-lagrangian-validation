#!/usr/bin/env python3
"""Draw the initial state of every Benchmark 2 test (no solver run).

Writes <test>/figures/<X>_setup.png: the domain (side view) with its mesh, the
particles at their start positions with their start velocity, and the air.

    python3 draw_setups.py          (from the benchmark folder)
"""
import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
COL = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]


def domain(ax, Lx, Ly, nx=1, ny=1, fc="#eef4fb", ec="#6c8ebf"):
    ax.add_patch(Rectangle((0, 0), Lx, Ly, fill=True, fc=fc, ec=ec, lw=1.5))
    for i in range(1, nx):
        ax.plot([i * Lx / nx] * 2, [0, Ly], color="#b8cbe4", lw=0.6)
    for j in range(1, ny):
        ax.plot([0, Lx], [j * Ly / ny] * 2, color="#b8cbe4", lw=0.6)


def arrow(ax, x, y, dx, dy, color, text=None, fs=9, **kw):
    ax.annotate("", xy=(x + dx, y + dy), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8))
    if text:
        ax.text(x + dx, y + dy, " " + text, color=color, fontsize=fs, va="center", **kw)


def curved(ax, x, y, dx, dy, color, rad, text=None, fs=9):
    ax.add_patch(FancyArrowPatch((x, y), (x + dx, y + dy), color=color, lw=1.4,
                                 arrowstyle="-|>", mutation_scale=14,
                                 connectionstyle=f"arc3,rad={rad}"))
    if text:
        ax.text(x + dx / 2, y + dy / 2, text, color=color, fontsize=fs,
                ha="center", va="center")


def particle(ax, x, y, color, label, dy_label):
    ax.plot(x, y, "o", color=color, ms=8, zorder=5)
    ax.text(x, y + dy_label, label, color=color, fontsize=8.5, ha="center",
            va="bottom" if dy_label > 0 else "top")


def save(fig, test, name):
    d = os.path.join(HERE, test, "figures")
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {test}/figures/{name}")


# ─── E: coupled uniform flow ──────────────────────────────────────────
def draw_E():
    fig, ax = plt.subplots(figsize=(8, 6))
    domain(ax, 1, 1, 10, 10)
    # inlet / outlet
    ax.plot([0, 0], [0, 1], color="#2e7d32", lw=3)
    ax.plot([1, 1], [0, 1], color="#2e7d32", lw=3)
    ax.text(-0.04, 0.5, "inlet\nU = (2 0 0)", color="#2e7d32", fontsize=9,
            ha="right", va="center")
    ax.text(1.04, 0.5, "outlet\np = 0", color="#2e7d32", fontsize=9,
            ha="left", va="center")
    for yy in (0.15, 0.85):
        arrow(ax, 0.05, yy, 0.18, 0, "#6c8ebf")
    ax.text(0.26, 0.85, "air: uniform wind U = (2, 0, 0) m/s\n(solved every step: PISO)",
            color="#4a6fa5", fontsize=9.5, va="center")

    ax.plot(0.2, 0.5, "o", color=COL[1], ms=11, zorder=4)
    ax.plot(0.2, 0.5, "o", color=COL[0], ms=6, zorder=5)
    ax.plot(0.9, 0.5, "o", color=COL[2], ms=8, zorder=5)
    ax.text(0.2, 0.55, "P0, P1", color="k", fontsize=8.5, ha="center", va="bottom")
    ax.text(0.9, 0.55, "P2", color=COL[2], fontsize=8.5, ha="center", va="bottom")
    ax.text(
        0.04, 0.34,
        "P0: d = 50 µm  u0 = (0, 0, 0)   x0 = (0.2, 0.5, 0.5)\n"
        "P1: d = 40 µm  u0 = (0, 1, 0)   x0 = (0.2, 0.5, 0.5)\n"
        "    (pushed sideways in y, out of the picture)\n"
        "P2: d = 60 µm  u0 = (0, 0, 0)   x0 = (0.9, 0.5, 0.5)",
        fontsize=7.8, va="top", family="monospace",
    )
    arrow(ax, 0.93, 0.8, 0, -0.1, "k", "g")
    ax.text(0.02, 0.02, "ρp = 2500 kg/m³, stokesDrag + gravity, dt = 5e-4 s",
            fontsize=9.5)
    ax.set(xlim=(-0.33, 1.33), ylim=(-0.03, 1.03), aspect="equal",
           xticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0], yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
           xlabel="x [m]", ylabel="y [m]")
    ax.set_title("Test 5 – initial state\nbox 1 × 1 × 1 m, 10 × 10 × 10 cells, "
                 "sides symmetry (uniform flow is exact)", fontsize=11)
    save(fig, "05_uniform_flow", "E_setup.png")


# ─── F: live non-uniform air (lid-driven cavity) ──────────────────────
def draw_F():
    fig, ax = plt.subplots(figsize=(7.5, 7))
    domain(ax, 1, 1, 20, 20)

    # moving lid
    ax.plot([0, 1], [1, 1], color="tab:red", lw=4, solid_capstyle="butt")
    arrow(ax, 0.38, 1.055, 0.24, 0, "tab:red")
    ax.text(0.5, 1.10, "moving lid: U = (1, 0, 0) m/s", color="tab:red",
            fontsize=10, ha="center", va="bottom")
    # fixed walls
    for x0, y0, x1, y1 in ((0, 0, 0, 1), (1, 0, 1, 1), (0, 0, 1, 0)):
        ax.plot([x0, x1], [y0, y1], color="k", lw=3, solid_capstyle="butt")
    ax.text(0.5, -0.07, "fixed walls: no-slip", fontsize=9.5, ha="center", va="top")
    ax.text(-0.03, 0.5, "fixed wall\n(no-slip)", fontsize=9, ha="right", va="center",
            rotation=90)

    # the air: at rest at t = 0, then swirls
    ax.text(0.5, 0.62, "air at rest at t = 0:  U = (0, 0, 0)\n"
                       "solved every step (PISO) → swirl",
            fontsize=10, ha="center", va="center",
            bbox=dict(fc="white", ec="#6c8ebf", lw=0.8, boxstyle="round,pad=0.4"))
    curved(ax, 0.30, 0.42, 0.40, 0.0, "#4a6fa5", -0.55)
    ax.text(0.50, 0.30, "U(t)", color="#4a6fa5", fontsize=10, ha="center")

    # particles near the lid
    xs = [0.113, 0.313, 0.513, 0.713, 0.913]
    for i, (x, c) in enumerate(zip(xs, COL)):
        particle(ax, x, 0.913, c, f"P{i}", -0.022)
    ax.text(0.5, 0.83, "5 particles, released at rest: u0 = (0, 0, 0)", fontsize=9.5,
            ha="center", va="top")

    ax.text(0.02, 0.02,
            "d = 1 mm, ρp = 2500 kg/m³, τp = 0.0116 s\n"
            "stokesDrag only (no gravity), RK2\n"
            "dt = 5e-3 s, 600 steps to t = 3 s\n"
            "z: one cell (0.1 m), front/back empty",
            fontsize=9, va="bottom",
            bbox=dict(fc="white", ec="#cccccc", lw=0.6, boxstyle="round,pad=0.35"))

    ax.set(xlim=(-0.12, 1.02), ylim=(-0.12, 1.16), aspect="equal",
           xlabel="x [m]", ylabel="y [m]")
    ax.set_title("Test 6 – initial state\ncavity 1 × 1 × 0.1 m, 20 × 20 × 1 cells "
                 "(cell 0.05 m), closed", fontsize=11)
    save(fig, "06_cavity", "F_setup.png")


if __name__ == "__main__":
    draw_E()
    draw_F()
