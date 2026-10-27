#!/usr/bin/env python3
"""Test B - forces and drag: settling of five solid particles from rest.

  results/stokes_<scheme>.csv   stokesDrag: compare u_z(t) with the exact solution (every diameter)
  results/sphere_<scheme>.csv   sphereDrag (Putnam):   final u_z vs the terminal velocity of the force balance
  results/SN_<scheme>.csv       SchillerNaumannDrag:   same            (<scheme> = RK2, Euler)

Prints the tables and PASS or FAIL, writes figures/B_terminal.png, B_settling.png and B_dragRe.png.
Run from the test folder:   python3 scripts/plot_B.py
"""
import os

import numpy as np
import matplotlib.pyplot as plt

# case parameters (must match the case files)
rho_p = 2500.0               # particle density [kg/m^3]
rho_f, nu = 1.2, 1.5e-5      # air
g = 9.81                     # |g| [m/s^2], gravity along -z
mu = rho_f * nu
g_eff = g * (1 - rho_f / rho_p)  # gravity with buoyancy


def CdRe_putnam(Re):
    return np.where(Re > 1000, 0.424 * Re, 24 * (1 + Re ** (2 / 3) / 6))


def CdRe_SN(Re):
    return np.where(Re > 1000, 0.44 * Re, 24 * (1 + 0.15 * Re**0.687))


def CdRe_stokes(Re):
    return 24.0 + 0 * Re


def tau_p(d):
    return rho_p * d**2 / (18 * mu)


def terminal_velocity(d, CdRe):
    """force balance (rho_p - rho_f) V g = (pi/8) mu d CdRe(Re) v, by fixed-point iteration"""
    v_stokes = tau_p(d) * g_eff
    v = v_stokes
    for _ in range(1000):
        Re = rho_f * v * d / mu
        v = v_stokes * 24 / CdRe(Re)
    return v


def read(name):
    f = f"results/{name}.csv"
    if not os.path.exists(f):
        return None
    D = np.genfromtxt(f, delimiter=",", names=True)
    return {dd: D[D["d"] == dd] for dd in np.unique(D["d"])}


os.makedirs("figures", exist_ok=True)
ok_all = True

SCHEMES = ["RK2", "Euler"]
STOKES_LIMIT = {"RK2": 1e-3, "Euler": 2e-2}   # Euler is 1st order: larger error while speeding up

# ─── 1. Stokes drag: exact transient for every diameter, both integrators ───
for scheme in SCHEMES:
    S = read(f"stokes_{scheme}")
    if not S:
        continue
    print(f"1. STOKES DRAG vs EXACT, {scheme}   u_z(t) = -v_s (1 - exp(-t/tau_p)),  v_s = tau_p g (1 - rho_f/rho_p)")
    print(f"{'d [um]':>7s} {'tau_p [s]':>10s} {'v_s [m/s]':>10s} {'max|u-exact|/v_s':>17s}  verdict (limit {STOKES_LIMIT[scheme]:g})")
    for dd, D in sorted(S.items()):
        tp = tau_p(dd)
        vs = tp * g_eff
        ue = -vs * (1 - np.exp(-D["t"] / tp))
        e = np.abs(D["uz"] - ue).max() / vs
        ok = e < STOKES_LIMIT[scheme]
        ok_all &= ok
        print(f"{dd * 1e6:7.0f} {tp:10.3e} {vs:10.4e} {e:17.3e}  {'PASS' if ok else 'FAIL'}")
    print()

# ─── 2. Nonlinear drag: terminal velocity, both integrators ───
models = [("sphere", "sphereDrag (Putnam)", CdRe_putnam), ("SN", "Schiller-Naumann", CdRe_SN)]
print("2. TERMINAL VELOCITY vs FORCE BALANCE  (u_z at t_end)")
print(f"{'drag':8s} {'scheme':6s} {'d [um]':>7s} {'v_t ref [m/s]':>14s} {'v_t NeoFOAM':>14s} {'rel. error':>11s} {'Re_t':>8s}  verdict")
fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
for a, (n, lbl, CdRe) in zip(ax, models):
    runs = {s: read(f"{n}_{s}") for s in SCHEMES}
    sizes = sorted(next(r for r in runs.values() if r))
    theory = [terminal_velocity(dd, CdRe) for dd in sizes]
    x = np.arange(len(sizes))
    w = 0.27
    bars = a.bar(x - w, theory, w, color="lightgray", edgecolor="k", lw=0.5, label="theory (gravity = drag)")
    for b, v in zip(bars, theory):
        a.text(b.get_x() + 1.5 * w, v * 1.06, f"{v:.3g} m/s", ha="center", fontsize=8)
    for k, (s, col) in enumerate(zip(SCHEMES, ["tab:blue", "tab:green"])):
        if not runs[s]:
            continue
        vals = []
        for dd, vt in zip(sizes, theory):
            vn = -runs[s][dd]["uz"][-1]
            e = abs(vn - vt) / vt
            ok = e < 1e-4
            ok_all &= ok
            print(f"{n:8s} {s:6s} {dd * 1e6:7.0f} {vt:14.6e} {vn:14.6e} {e:11.2e} {rho_f * vt * dd / mu:8.3f}  {'PASS' if ok else 'FAIL'}")
            vals.append(vn)
        a.bar(x + k * w, vals, w, color=col, edgecolor="k", lw=0.5, label=f"NeoFOAM {s}")
    a.set_yscale("log")
    a.set_xticks(x, [f"{dd * 1e6:.0f} µm" for dd in sizes])
    a.set(xlabel="particle diameter", title=lbl, ylim=(4e-3, 8))
    a.legend(fontsize=8, loc="upper left")
ax[0].set_ylabel("final falling speed [m/s]  (log scale)")
fig.suptitle("Test B – final falling speed of each size: theory, NeoFOAM RK2 and NeoFOAM Euler are the same")
fig.tight_layout()
fig.savefig("figures/B_terminal.png", dpi=150)

# ─── 3. Settling with sphereDrag: falling speed over time ───
P = read("sphere_RK2")
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
if P:
    def reference(d, n_steps=4000):
        """fine RK4 of du/dt = g_eff - drag/m (Putnam), from rest; returns t, u (downwards)"""
        vt = terminal_velocity(d, CdRe_putnam)
        t_end = min(8.0, 15 * vt / g_eff)          # long enough to settle
        dt = t_end / n_steps
        m = rho_p * np.pi * d**3 / 6

        def a(u):
            Re = rho_f * d * abs(u) / mu
            return g_eff - np.pi / 8 * mu * d * float(CdRe_putnam(Re)) * u / m

        ts, us, u = [0.0], [0.0], 0.0
        for k in range(n_steps):
            k1 = a(u); k2 = a(u + dt / 2 * k1); k3 = a(u + dt / 2 * k2); k4 = a(u + dt * k3)
            u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            ts.append((k + 1) * dt); us.append(u)
        return np.array(ts), np.array(us), vt

    fig, ax = plt.subplots(figsize=(8, 5))
    for (dd, D), c in zip(sorted(P.items()), colors):
        tr, ur, vt = reference(dd)
        ax.plot(D["t"], -D["uz"], "-", color=c, lw=2, label=f"d = {dd * 1e6:.0f} µm")
        ax.plot(tr, ur, "k--", lw=0.8)
        ax.text(2.55, vt, f"{vt:.3g} m/s", color=c, fontsize=8, va="center")
    ax.plot([], [], "k--", lw=0.8, label="theory (fine reference)")
    ax.set(xlabel="t [s]", ylabel="falling speed [m/s]", xlim=(0, 2.5),
           title="Falling from rest until drag = gravity (sphereDrag, RK2)")
    ax.legend(fontsize=8, loc="center right")
    fig.suptitle("Test B – five solid particles falling from rest (ρp = 2500 kg/m³)")
    fig.tight_layout()
    fig.savefig("figures/B_settling.png", dpi=150)

# ─── 4. Drag law vs Reynolds number, from the NeoFOAM terminal velocities ───
# at the terminal velocity drag = gravity, so (Cd Re / 24) = v_stokes / v_t: a drag factor
# measured from NeoFOAM, to compare with the curve of the drag law at Re_t = rho_f v_t d / mu
fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
Re_line = np.logspace(-3, np.log10(300), 400)
print("\n4. DRAG FACTOR (Cd Re / 24) FROM THE NeoFOAM TERMINAL VELOCITY")
print(f"{'drag':8s} {'scheme':6s} {'d [um]':>7s} {'Re_t':>8s} {'factor NeoFOAM':>15s} {'factor law':>11s}")
for a, (n, lbl, CdRe) in zip(ax, models):
    a.plot(Re_line, CdRe(Re_line) / 24, "-", color="tab:gray", lw=2, label=f"{lbl} law")
    a.plot(Re_line, CdRe_stokes(Re_line) / 24, "k:", lw=1, label="Stokes (= 1)")
    for s, style in [("RK2", "o"), ("Euler", "x")]:
        runs = read(f"{n}_{s}")
        if not runs:
            continue
        for (dd, D), c in zip(sorted(runs.items()), colors):
            vt = -D["uz"][-1]
            Re_t = rho_f * vt * dd / mu
            factor = tau_p(dd) * g_eff / vt
            print(f"{n:8s} {s:6s} {dd * 1e6:7.0f} {Re_t:8.3f} {factor:15.4f} {float(CdRe(Re_t)) / 24:11.4f}")
            a.plot(Re_t, factor, style, color=c, ms=10 if s == "RK2" else 8, mfc="none" if s == "RK2" else None,
                   label=f"d = {dd * 1e6:.0f} µm" if s == "RK2" else None)
    a.plot([], [], "ko", ms=10, mfc="none", label="NeoFOAM RK2")
    a.plot([], [], "kx", ms=8, label="NeoFOAM Euler")
    a.set_xscale("log")
    a.set(xlabel="particle Reynolds number Re at the final speed", title=lbl)
    a.legend(fontsize=8, loc="upper left")
ax[0].set_ylabel("drag compared with Stokes drag  (Cd Re / 24)")
fig.suptitle("Test B – each particle feels the drag of its Reynolds number: NeoFOAM points on the drag law")
fig.tight_layout()
fig.savefig("figures/B_dragRe.png", dpi=150)

print(f"\nOVERALL: {'PASS' if ok_all else 'FAIL'}   (figures in figures/)")
