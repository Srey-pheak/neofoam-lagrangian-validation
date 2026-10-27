#!/usr/bin/env python3
"""Test C - 3D motion and boundaries: three particles in a uniform wind with gravity, Stokes drag.

Exact solution per component (Stokes drag is linear):
    u_inf = U_f + tau_p g (1 - rho_f/rho_p)
    u(t)  = u_inf + (u0 - u_inf) exp(-t/tau_p)
    x(t)  = x0 + u_inf t + (u0 - u_inf) tau_p (1 - exp(-t/tau_p))

  particles inside the box: max error of u and x (all components) vs exact
  particle leaving through +x (x = 1 m): exit time vs the exact crossing time, and the
  particle must stop moving once it has left (deactivated)

Both integrators: results/C_RK2.csv and results/C_Euler.csv (limit 1e-3 for RK2, 2e-2 for Euler,
which is 1st order). Prints the results and PASS or FAIL, writes figures/C_trajectories_<scheme>.png
(paths from the side and from above, and u_x, u_y, u_z over time: NeoFOAM line, exact dashed).
Run from the test folder:   python3 scripts/plot_C.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SCHEMES = [s for s in ("RK2", "Euler") if os.path.exists(f"results/C_{s}.csv")]
if not SCHEMES:
    sys.exit("need results/C_RK2.csv, results/C_Euler.csv: run scripts/run_C.sh")
LIMIT = {"RK2": 1e-3, "Euler": 2e-2}   # Euler is 1st order: larger error
COLORS = ["tab:blue", "tab:orange", "tab:green"]

# case parameters (must match the case files)
rho_p = 2500.0
rho_f, nu = 1.2, 1.5e-5
mu = rho_f * nu
Uf = np.array([2.0, 0.0, 0.0])          # wind [m/s]
g = np.array([0.0, 0.0, -9.81])
L = 1.0                                 # box size [m], the exit face is x = L
particles = {                           # id (injection order): d, x0, u0
    0: (50e-6, [0.2, 0.5, 0.5], [0.0, 0.0, 0.0]),
    1: (40e-6, [0.2, 0.5, 0.5], [0.0, 1.0, 0.0]),
    2: (60e-6, [0.9, 0.5, 0.5], [0.0, 0.0, 0.0]),
}


def exact(t, d, x0, u0):
    tau = rho_p * d**2 / (18 * mu)
    uinf = Uf + tau * g * (1 - rho_f / rho_p)
    e = np.exp(-t / tau)[:, None]
    u = uinf + (np.array(u0) - uinf) * e
    x = np.array(x0) + uinf * t[:, None] + (np.array(u0) - uinf) * tau * (1 - e)
    return x, u, tau


os.makedirs("figures", exist_ok=True)
ok_all = True

for scheme in SCHEMES:
    print(f"=== {scheme} (limit {LIMIT[scheme]:g}) ===")
    D = np.genfromtxt(f"results/C_{scheme}.csv", delimiter=",", names=True)
    dt = np.diff(np.unique(D["t"]))[0]
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 6)
    ax = [fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])]
    av = [fig.add_subplot(gs[1, 2 * k:2 * k + 2]) for k in range(3)]

    for pid, (d, x0, u0) in particles.items():
        P = D[D["id"] == pid]
        if len(P) == 0 or abs(P["d"][0] - d) > 1e-12:
            print(f"particle {pid}: not found or wrong diameter in the CSV")
            ok_all = False
            continue
        t = P["t"]
        X = np.column_stack([P["x"], P["y"], P["z"]])
        U = np.column_stack([P["ux"], P["uy"], P["uz"]])
        xe, ue, tau = exact(t, d, x0, u0)
        inside = xe[:, 0] < L                       # exact solution still inside the box
        uscale = np.abs(ue).max()

        mask = inside                               # times to plot the error at
        if inside.all():
            e_u = np.abs(U - ue).max() / uscale
            e_x = np.abs(X - xe).max()
            ok = e_u < LIMIT[scheme]
            print(f"particle {pid} (d = {d * 1e6:.0f} µm, tau_p = {tau:.4f} s, u0 = {u0}):")
            print(f"    max |u - exact| / max|u| = {e_u:.2e},  max |x - exact| = {e_x:.2e} m   ->  {'PASS' if ok else 'FAIL'}")
        else:
            # exact crossing time of x = L (linear interpolation on a fine grid)
            tf = np.linspace(0, t[-1], 200001)
            xf, _, _ = exact(tf, d, x0, u0)
            t_exit = tf[np.argmax(xf[:, 0] >= L)]
            # NeoFOAM: first output time with x >= L, and the particle must not move afterwards
            k = np.argmax(X[:, 0] >= L)
            t_out = t[k]
            frozen = np.abs(X[k:] - X[k]).max() == 0.0 and np.abs(U[k:] - U[k]).max() == 0.0
            before = t < t_exit - dt
            mask = before                           # after leaving it is frozen: no exact solution to compare
            e_u = np.abs(U[before] - ue[before]).max() / uscale
            ok = frozen and 0 <= t_out - t_exit <= 2 * dt and e_u < LIMIT[scheme]
            print(f"particle {pid} (d = {d * 1e6:.0f} µm, tau_p = {tau:.4f} s) leaves through x = {L} m:")
            print(f"    exact crossing t = {t_exit:.5f} s, NeoFOAM first x >= {L}: t = {t_out:.5f} s "
                  f"(+{(t_out - t_exit) / dt:.1f} dt), frozen afterwards: {frozen}")
            print(f"    before leaving: max |u - exact| / max|u| = {e_u:.2e}   ->  {'PASS' if ok else 'FAIL'}")
        ok_all &= ok

        lbl = f"P{pid}, d = {d * 1e6:.0f} µm"
        c = COLORS[pid]
        for a, (i, j) in zip(ax, [(0, 2), (0, 1)]):
            a.plot(X[:, i], X[:, j], "-", color=c, lw=2.5, label=lbl)
            a.plot(xe[inside, i], xe[inside, j], "k--", lw=0.8)
        for k, a in enumerate(av):                  # until it leaves: afterwards it is frozen
            a.plot(t[mask], U[mask, k], "-", color=c, lw=2.5, label=lbl)
            a.plot(t[mask], ue[mask, k], "k--", lw=0.8)

    for a in ax:
        a.axvline(L, color="k", lw=1)
        a.text(L, a.get_ylim()[0], " exit face x = 1 m", rotation=90, va="bottom", fontsize=8)
    ax[0].set(xlabel="x [m]", ylabel="z [m]", title="Side view (x–z): wind →, gravity ↓")
    ax[1].set(xlabel="x [m]", ylabel="y [m]", title="Top view (x–y)")
    ax[0].plot([], [], "k--", lw=0.8, label="exact")
    ax[0].legend(fontsize=8)
    for a, n in zip(av, "xyz"):
        a.set(xlabel="t [s]", ylabel=f"u_{n} [m/s]", title=f"Velocity in {n}")
    av[0].plot([], [], "k--", lw=0.8, label="exact")
    av[0].legend(fontsize=8)
    fig.suptitle(f"Test C – {scheme}: uniform wind (2 0 0) m/s + gravity, Stokes drag, 10×10×10 cells")
    fig.tight_layout()
    fig.savefig(f"figures/C_trajectories_{scheme}.png", dpi=150)
    print()

print(f"OVERALL: {'PASS' if ok_all else 'FAIL'}   (figures/C_trajectories_<scheme>.png)")
