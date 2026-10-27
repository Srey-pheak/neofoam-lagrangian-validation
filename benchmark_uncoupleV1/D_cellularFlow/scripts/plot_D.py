#!/usr/bin/env python3
"""Test D - cellular flow (Maxey 1987): 4 particle sizes (St = 0.1, 0.5, 1, 2) in 2 x 2 vortices.

Two references, both computed here independently of NeoFOAM:

  1. same-cell reference: the same steps as the solver (each step: find the cell, take the
     cell value of U, one RK2 step with Stokes drag), but the cell is found with arithmetic,
     cell = floor(x / dx) + nx floor(y / dx), instead of OpenFOAM's findCell.
     If NeoFOAM finds the correct cell in every step, the two agree to round-off.
     PASS: max |x_NeoFOAM - x_ref| < 1e-9 m (the cell size is 3e-4 m) and the same
     particles leave the domain.
  2. smooth-flow reference: fine RK4 with the exact (smooth) cellular flow. It shows how much
     the cell value of U (a step-like field) changes the paths. Information only, no pass/fail.

Prints the tables and PASS or FAIL, writes figures/D_particles.png (positions at t = 0, 0.1, 0.2 s).
Run from the test folder:   python3 scripts/plot_D.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from params import *  # noqa: E402,F403

if not (os.path.exists("results/D.csv") and os.path.exists("results/U_cells")):
    sys.exit("need results/D.csv and results/U_cells: run scripts/run_D.sh")

D = np.genfromtxt("results/D.csv", delimiter=",", names=True)
times = np.unique(D["t"])
ids = np.unique(D["id"]).astype(int)
n = len(ids)
# NeoFOAM data as arrays [output time, particle]
X = np.full((len(times), n, 3), np.nan)
for k, c in enumerate("xyz"):
    X[:, :, k] = D[c].reshape(len(times), n)
d = D["d"][:n]


def read_cell_values(path="results/U_cells"):
    """the cell values of U that NeoFOAM used (0/U of the run, written by setExprFields)"""
    lines = open(path).read().split("internalField")[1].split("(", 1)[1].split("\n)")[0]
    vals = [list(map(float, s.strip("() ").split())) for s in lines.strip().splitlines() if s.strip().startswith("(")]
    return np.array(vals)


Ucell = read_cell_values()
x0 = X[0].copy()
mu_f = rho_f * nu
m = rho_p * np.pi / 6.0 * d * d * d
k_drag = np.pi / 8.0 * mu_f * d * 24.0             # Stokes: F = k (uf - u)
n_steps = int(round(end_time / dt))


def inside(x):
    return (x[:, 0] >= 0) & (x[:, 0] < 2 * L) & (x[:, 1] >= 0) & (x[:, 1] < 2 * L) & (x[:, 2] >= 0) & (x[:, 2] < dz)


# ─── 1. same-cell reference: the solver's steps, cell found by arithmetic ───
x = x0.copy()
u = np.zeros_like(x)
active = np.ones(n, bool)
ref = [x.copy()]
for step in range(1, n_steps + 1):
    # locate (Q1, Q2, Q10)
    active &= inside(x)
    i = np.clip(np.floor(x[:, 0] / dx).astype(int), 0, nx - 1)
    j = np.clip(np.floor(x[:, 1] / dx).astype(int), 0, ny - 1)
    uf = Ucell[i + nx * j]
    # evolve (Q3 - Q8): RK2 with the cell value of U, Stokes drag
    a = lambda up: (uf - up) * k_drag[:, None] / m[:, None]   # noqa: E731
    a1 = a(u)
    u_star = u + a1 * dt
    x_star = x + u * dt
    a2 = a(u_star)
    u_new = u + (a1 + a2) * (0.5 * dt)
    x_new = x + (u + u_star) * (0.5 * dt)
    u = np.where(active[:, None], u_new, u)
    x = np.where(active[:, None], x_new, x)
    if step % write_interval == 0:
        ref.append(x.copy())
ref = np.array(ref)
left_ref = ~active

# ─── 2. smooth-flow reference: fine RK4 with the exact flow ───
def flow_v(p):
    return np.stack([U0 * np.sin(np.pi * p[:, 0] / L) * np.cos(np.pi * p[:, 1] / L),
                     -U0 * np.cos(np.pi * p[:, 0] / L) * np.sin(np.pi * p[:, 1] / L),
                     0 * p[:, 0]], axis=1)


tau_p = rho_p * d**2 / (18 * mu_f)
sub = 4                                            # 4 RK4 steps per solver step
h = dt / sub
xs, us = x0.copy(), np.zeros_like(x0)
act = np.ones(n, bool)
smooth = [xs.copy()]
for step in range(1, n_steps * sub + 1):
    act &= inside(xs)
    f = lambda p, v: ((flow_v(p) - v) / tau_p[:, None], v)    # noqa: E731
    k1u, k1x = f(xs, us)
    k2u, k2x = f(xs + h / 2 * k1x, us + h / 2 * k1u)
    k3u, k3x = f(xs + h / 2 * k2x, us + h / 2 * k2u)
    k4u, k4x = f(xs + h * k3x, us + h * k3u)
    un = us + h / 6 * (k1u + 2 * k2u + 2 * k3u + k4u)
    xn = xs + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
    us = np.where(act[:, None], un, us)
    xs = np.where(act[:, None], xn, xs)
    if step % (write_interval * sub) == 0:
        smooth.append(xs.copy())
smooth = np.array(smooth)

# NeoFOAM: a particle that left stays frozen at its last position, which is outside the domain
left_neo = ~inside(X[-1])

# ─── table ───
ok_all = True
print("TEST D - cellular flow, 2 x 2 vortices, "
      f"{nx} x {ny} cells, U0 = {U0:g} m/s, L = {L:g} m, tau_f = {tau_f:g} s, {n_steps} steps\n")
print(f"{'St':>4s} {'d [um]':>7s} {'n':>4s} {'left NeoFOAM':>13s} {'left ref':>9s} "
      f"{'max|x - same-cell ref| [m]':>27s} {'max|x - smooth ref| [m]':>24s}  verdict")
groups = []
for St in ST:
    g = np.where(np.isclose(d, diameter(St), rtol=1e-12))[0]
    groups.append((St, g))
    dx_same = np.nanmax(np.abs(X[:, g] - ref[:, g]))
    dx_smooth = np.nanmax(np.abs(X[:, g] - smooth[:, g]))
    same_left = np.array_equal(left_neo[g], left_ref[g])
    ok = dx_same < 1e-9 and same_left
    ok_all &= ok
    print(f"{St:4g} {diameter(St) * 1e6:7.1f} {len(g):4d} {left_neo[g].sum():13d} {left_ref[g].sum():9d} "
          f"{dx_same:27.2e} {dx_smooth:24.2e}  {'PASS' if ok else 'FAIL'}")
print(f"\n(cell size dx = {dx:.2e} m; the smooth-flow column is information only: "
      "the solver uses the cell value of U)")

# ─── how many particles are still inside the domain ───
T_SHOW = [0.1, 0.2, 0.3, 0.5, 1.0]
k_show = [int(round(t / (dt * write_interval))) for t in T_SHOW]
print("\nPARTICLES STILL INSIDE THE DOMAIN (of 100), NeoFOAM | smooth-flow reference")
print(f"{'St':>4s}  " + "  ".join(f"{'t = ' + format(t, 'g') + ' s':>11s}" for t in T_SHOW))
for St, g in groups:
    print(f"{St:4g}  " + "  ".join(f"{int(inside(X[k, g]).sum()):4d} | {int(inside(smooth[k, g]).sum()):4d}"
                                    for k in k_show))

# ─── figure: positions at t = 0, 1 and 2 tau_f, one row per St, on the streamlines ───
T_FIG = [0.0, 0.1, 0.2]
gx = np.linspace(0, 2 * L, 200)
GX, GY = np.meshgrid(gx, gx)
psi = U0 * L / np.pi * np.sin(np.pi * GX / L) * np.sin(np.pi * GY / L)   # stream function
fig, ax = plt.subplots(len(ST), len(T_FIG), figsize=(10, 13.5))
for r, ((St, g), c) in enumerate(zip(groups, ["tab:blue", "tab:orange", "tab:green", "tab:red"])):
    for col, t in enumerate(T_FIG):
        a = ax[r, col]
        k = int(round(t / (dt * write_interval)))
        a.contour(GX * 1e3, GY * 1e3, psi, levels=12, colors="lightgray", linewidths=0.8)
        inn = inside(X[k, g])
        a.plot(X[k, g][inn, 0] * 1e3, X[k, g][inn, 1] * 1e3, "o", color=c, ms=4)
        a.set(xlim=(0, 20), ylim=(0, 20), aspect="equal", xticks=[0, 10, 20], yticks=[0, 10, 20])
        a.set_title(f"St = {St:g} (d = {diameter(St) * 1e6:.0f} µm), t = {t:g} s: {inn.sum()} inside", fontsize=9)
        if r == len(ST) - 1:
            a.set_xlabel("x [mm]")
        if col == 0:
            a.set_ylabel("y [mm]")
fig.suptitle(f"Test D – cellular flow (4 vortices, τf = {tau_f:g} s): particles are thrown out of the vortices\n"
             "grey lines: streamlines; a particle that reaches the domain edge leaves")
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/D_particles.png", dpi=150)
print(f"\nOVERALL: {'PASS' if ok_all else 'FAIL'}   (figures/D_particles.png)")
