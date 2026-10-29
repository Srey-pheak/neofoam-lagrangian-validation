#!/usr/bin/env python3
"""Test 1 - integrators: NeoFOAM particle velocity vs the exact Stokes-drag solution.

  order of accuracy: results/<scheme>_dt<dt>.csv  -> figures/A_velocity.png, A_convergence.png

Prints the error / order table and PASS or FAIL. Run from the test folder:
    python3 scripts/plot_A.py
"""
import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt

# case parameters (must match the case files)
rho_p, d = 2500.0, 5e-5      # particle: density [kg/m^3], diameter [m]
rho_f, nu = 1.2, 1.5e-5      # air: density [kg/m^3], kinematic viscosity [m^2/s]
u0 = 0.1                     # initial velocity [m/s]
mu = rho_f * nu
tau = rho_p * d**2 / (18 * mu)


def exact_u(t):
    return u0 * np.exp(-t / tau)


def exact_dx(t):
    """x(t) - x0"""
    return u0 * tau * (1 - np.exp(-t / tau))


def read_runs(pattern):
    """{scheme: [(dt, data), ...]} from files <scheme>_dt<dt>.csv, largest dt first"""
    runs = {}
    for f in glob.glob(pattern):
        m = re.match(r"(\w+?)_dt([0-9.eE+-]+)\.csv", os.path.basename(f))
        if not m:
            continue
        D = np.genfromtxt(f, delimiter=",", names=True)
        D = D[D["id"] == 0]
        runs.setdefault(m.group(1), []).append((float(m.group(2)), D))
    for lst in runs.values():
        lst.sort(key=lambda r: -r[0])
    return runs


os.makedirs("figures", exist_ok=True)
markers = {"Euler": "o", "RK2": "s"}
expected = {"Euler": 1, "RK2": 2}

# ─────────────────────────────────────────────────────────────────────
# 1. Order of accuracy
# ─────────────────────────────────────────────────────────────────────
runs = read_runs("results/*_dt*.csv")
print(f"tau_p = {tau:.6f} s\n")
print("1. ORDER OF ACCURACY (error = max over t, vs exact)")
print(f"{'scheme':6s} {'dt [s]':>10s} {'dt/tau':>7s} {'e_u':>11s} {'p_u':>5s} {'e_x':>11s} {'p_x':>5s}")
conv = {}
for scheme, lst in sorted(runs.items()):
    conv[scheme] = []
    prev = None
    for dt, D in lst:
        t = D["t"]
        x0 = D["x"][0]
        e_u = np.abs(D["ux"] - exact_u(t)).max() / u0
        e_x = np.abs(D["x"] - x0 - exact_dx(t)).max() / (u0 * tau)
        p_u = p_x = ""
        if prev:
            r = prev[0] / dt
            p_u = f"{np.log(prev[1] / e_u) / np.log(r):.2f}"
            p_x = f"{np.log(prev[2] / e_x) / np.log(r):.2f}"
        print(f"{scheme:6s} {dt:10.3e} {dt / tau:7.3f} {e_u:11.3e} {p_u:>5s} {e_x:11.3e} {p_x:>5s}")
        conv[scheme].append((dt, e_u, e_x, p_u))
        prev = (dt, e_u, e_x)
print("e_u = max|u - u_exact|/u0,  e_x = max|x - x_exact|/(u0 tau_p),  p = order between two dt")

# PASS: order between the two finest dt within 0.1 of the expected order
for scheme, c in sorted(conv.items()):
    if len(c) >= 2 and c[-1][3]:
        p = float(c[-1][3])
        ok = abs(p - expected.get(scheme, p)) < 0.1
        print(f"  {scheme}: finest order p_u = {p:.2f}, expected {expected.get(scheme)}  ->  {'PASS' if ok else 'FAIL'}")

if runs:
    # velocity over time: slows down to 0
    fig, ax = plt.subplots(figsize=(7, 4.5))
    tt = np.linspace(0, 0.2, 400)
    ax.plot(tt, exact_u(tt), "k-", lw=2, label="exact  $u_0 e^{-t/\\tau_p}$")
    for scheme, lst in sorted(runs.items()):
        for dt, D in lst:
            ax.plot(D["t"], D["ux"], "--", lw=1, label=f"{scheme}, dt = {dt:g} s")
    ax.axvline(tau, ls=":", c="gray")
    ax.text(tau, 0.09, "  t = τp: speed down to 37 %", color="gray", fontsize=8)
    ax.set(xlabel="t [s]", ylabel="$u_x$ [m/s]", title="Particle velocity", xlim=(0, 0.2))
    ax.legend(fontsize=7)
    fig.suptitle(f"Test 1 – Stokes drag relaxation (d = 50 µm, ρp = 2500 kg/m³, τp = {tau:.4f} s)")
    fig.tight_layout()
    fig.savefig("figures/A_velocity.png", dpi=150)

    # convergence
    fig, ax = plt.subplots(figsize=(6, 5))
    for scheme, c in sorted(conv.items()):
        dts = np.array([r[0] for r in c])
        eu = np.array([r[1] for r in c])
        ax.loglog(dts, eu, markers.get(scheme, "^") + "-", ms=7, label=scheme)
        p = expected.get(scheme)
        if p:
            ax.loglog(dts, eu[-1] * (dts / dts[-1]) ** p, "k:" if p == 1 else "k--", lw=1, label=f"slope {p}")
    ax.set(xlabel="dt [s]", ylabel="error = max $|u - u_{exact}|\\ /\\ u_0$", title="Test 1 – order of accuracy")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/A_convergence.png", dpi=150)

# ─────────────────────────────────────────────────────────────────────
# 2. Size of the time step: small, medium, big dt, normal axes
# ─────────────────────────────────────────────────────────────────────
steps = read_runs("results/timestep/*_dt*.csv")
if steps:
    print("\n2. SIZE OF THE TIME STEP")
    print(f"{'scheme':6s} {'dt [s]':>8s} {'dt/tau':>7s} {'u at t = 0.2 s':>15s}  verdict")
    colors = {0.001: "tab:blue", 0.02: "tab:orange", 0.05: "tab:red"}
    labels = {0.001: "small", 0.02: "medium", 0.05: "big"}
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    tt = np.linspace(0, 0.2, 400)
    for a, scheme in zip(ax, ["Euler", "RK2"]):
        a.plot(tt, exact_u(tt), "k-", lw=4, alpha=0.35, label="exact")
        for dt, D in sorted(steps.get(scheme, [])):
            ok = abs(D["ux"][-1]) < u0
            print(f"{scheme:6s} {dt:8g} {dt / tau:7.2f} {D['ux'][-1]:15.3e}  {'slows down (ok)' if ok else 'speeds up (WRONG)'}")
            style = "--" if dt < 0.005 else "o-"       # many points for the small dt: no markers
            a.plot(D["t"], D["ux"], style, ms=5, lw=1.5, color=colors.get(dt),
                   label=f"{labels.get(dt, '')} dt = {dt:g} s ({dt / tau:.2g} τp)")
        a.axhline(0, color="gray", lw=0.8)
        a.set(xlabel="t [s]", title=scheme, xlim=(0, 0.2), ylim=(-0.02, 0.2))
        a.legend()
    ax[0].set_ylabel("particle velocity u [m/s]")
    fig.suptitle("Test 1 – what a bigger time step does: Euler always slows down, RK2 speeds up when dt is too big",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig("figures/A_timestep.png", dpi=150)

print("\nfigures written to figures/")
