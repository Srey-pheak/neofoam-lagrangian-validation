#!/usr/bin/env python3
"""Draw the initial state of every Benchmark 1 test (no solver run).

Writes <test>/figures/<X>_setup.png: the domain (side view x-z) with its mesh, the particles at
their start positions with their start velocity, the air velocity and gravity.

    python3 draw_setups.py          (from the benchmark folder)
"""
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
COL = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]


def domain(ax, Lx, Lz, nx=1, nz=1):
    ax.add_patch(Rectangle((0, 0), Lx, Lz, fill=True, fc="#eef4fb", ec="#6c8ebf", lw=1.5))
    for i in range(1, nx):
        ax.plot([i * Lx / nx] * 2, [0, Lz], color="#b8cbe4", lw=0.6)
    for j in range(1, nz):
        ax.plot([0, Lx], [j * Lz / nz] * 2, color="#b8cbe4", lw=0.6)


def arrow(ax, x, z, dx, dz, color, text=None, **kw):
    ax.annotate("", xy=(x + dx, z + dz), xytext=(x, z),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8))
    if text:
        ax.text(x + dx, z + dz, " " + text, color=color, fontsize=9, va="center", **kw)


def particle(ax, x, z, color, label, dz_label):
    ax.plot(x, z, "o", color=color, ms=8, zorder=5)
    ax.text(x, z + dz_label, label, color=color, fontsize=8.5, ha="center",
            va="bottom" if dz_label > 0 else "top")


def save(fig, test, name):
    d = os.path.join(HERE, test, "figures")
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {test}/figures/{name}")


# ─── A: integrators ───────────────────────────────────────────────────
def draw_A():
    fig, ax = plt.subplots(figsize=(6.5, 6))
    L = 0.1
    domain(ax, L, L)
    particle(ax, 0.05, 0.05, "k", "d = 50 µm, ρp = 2500 kg/m³\nx0 = (0.05, 0.05, 0.05) m", -0.006)
    arrow(ax, 0.05, 0.05, 0.02, 0, "tab:red", "u0 = (0.1, 0, 0) m/s")
    ax.text(0.004, 0.095, "air at rest: U = (0, 0, 0)\nno gravity", fontsize=9.5, va="top")
    ax.set(xlim=(-0.005, 0.105), ylim=(-0.005, 0.105), aspect="equal", xlabel="x [m]", ylabel="z [m]")
    ax.set_title("Test 1 – initial state\nbox 0.1 × 0.1 × 0.1 m, 1 cell, all faces patch 'sides'", fontsize=11)
    save(fig, "01_integrators", "A_setup.png")


# ─── B: forces and drag ───────────────────────────────────────────────
def draw_B():
    fig, ax = plt.subplots(figsize=(6.5, 8))
    Lx, Lz, z0 = 0.5, 35.0, 34.5
    domain(ax, Lx, Lz)
    for d, x, c, s in zip([10, 50, 100, 200, 500], [0.05, 0.15, 0.25, 0.35, 0.45], COL, [1, -1, 1, -1, 1]):
        particle(ax, x, z0, c, f"{d} µm\nx = {x}", -1.0)
    ax.text(0.25, 30.5, "all released at rest: u0 = (0, 0, 0)\nstart y = 0.25 m, z = 34.5 m\nρp = 2500 kg/m³",
            ha="center", fontsize=9, va="top")
    arrow(ax, 0.25, 22, 0, -5, "k", "g = (0, 0, −9.81) m/s²")
    ax.text(0.02, 1.0, "air at rest: U = (0, 0, 0)", fontsize=9.5)
    ax.set(xlim=(-0.02, 0.52), ylim=(-0.8, 35.8), xlabel="x [m]  (stretched: not to scale)", ylabel="z [m]")
    ax.set_title("Test 2 – initial state\ncolumn 0.5 × 0.5 × 35 m, 1 cell, all faces patch 'sides'", fontsize=11)
    save(fig, "02_drag_settling", "B_setup.png")


# ─── C: 3D and boundaries ─────────────────────────────────────────────
def draw_C():
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    domain(ax, 1, 1, 10, 10)
    ax.plot([1, 1], [0, 1], color="tab:red", lw=3)
    ax.text(1.015, 0.5, "x = 1 m face", color="tab:red", fontsize=9, rotation=90, va="center")
    for zz in (0.12, 0.88):
        arrow(ax, 0.05, zz, 0.15, 0, "#6c8ebf")
    ax.text(0.22, 0.88, "air: uniform wind U = (2, 0, 0) m/s", color="#4a6fa5", fontsize=9.5, va="center")
    # P0 and P1 start at the same point; the labels are placed apart
    ax.plot(0.2, 0.5, "o", color=COL[1], ms=11, zorder=4)
    ax.plot(0.2, 0.5, "o", color=COL[0], ms=6, zorder=5)
    ax.text(0.22, 0.44, "P0: d = 50 µm, u0 = (0, 0, 0)", color=COL[0], fontsize=8.5, va="top")
    ax.text(0.22, 0.40, "P1: d = 40 µm, u0 = (0, 1, 0): sideways (y, out of the picture)", color=COL[1], fontsize=8.5, va="top")
    ax.text(0.22, 0.36, "both at x0 = (0.2, 0.5, 0.5)", color="k", fontsize=8.5, va="top")
    ax.plot(0.9, 0.5, "o", color=COL[2], ms=8, zorder=5)
    ax.text(0.88, 0.56, "P2: d = 60 µm, u0 = (0, 0, 0)\nx0 = (0.9, 0.5, 0.5)", color=COL[2], fontsize=8.5, ha="right", va="bottom")
    arrow(ax, 0.93, 0.8, 0, -0.1, "k", "g")
    ax.text(0.02, 0.02, "ρp = 2500 kg/m³", fontsize=9.5)
    ax.set(xlim=(-0.02, 1.08), ylim=(-0.02, 1.02), aspect="equal", xlabel="x [m]", ylabel="z [m]")
    ax.set_title("Test 3 – initial state\nbox 1 × 1 × 1 m, 10 × 10 × 10 cells, all faces patch 'sides'", fontsize=11)
    save(fig, "03_motion_3d_exit", "C_setup.png")


if __name__ == "__main__":
    draw_A()
    draw_B()
    draw_C()
