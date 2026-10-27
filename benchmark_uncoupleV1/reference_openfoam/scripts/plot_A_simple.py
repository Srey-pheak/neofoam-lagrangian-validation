#!/usr/bin/env python3
"""Case A (sphereDrag + Euler): velocity, position and order of accuracy, NeoFOAM and OpenFOAM.

Writes figures/ref_A_simple.png. Run from the reference_openfoam folder:
    python3 scripts/plot_A_simple.py
"""
import math

import numpy as np
import matplotlib.pyplot as plt

rho_p, rho_f, nu = 2500.0, 1.2, 1.5e-5
mu = rho_f * nu
d, u0, x0 = 5e-5, 0.1, 0.05
DT_SHOW = "1e-3"                                   # the run shown in the u(t) and x(t) plots
DTS = ["2e-3", "1e-3", "5e-4", "2.5e-4", "1.25e-4"]


def reference(t_end=0.2, dt=2e-6):
    """fine RK4 solution of the same equation (Putnam drag, air at rest)"""
    m = rho_p * math.pi * d**3 / 6

    def a(u):
        Re = rho_f * d * abs(u) / mu
        return -math.pi / 8 * mu * d * 24 * (1 + Re ** (2 / 3) / 6) * u / m

    ts, us, xs, u, x = [0.0], [u0], [x0], u0, x0
    for k in range(int(round(t_end / dt))):
        k1 = a(u); k2 = a(u + dt / 2 * k1); k3 = a(u + dt / 2 * k2); k4 = a(u + dt * k3)
        x += dt / 6 * (u + 2 * (u + dt / 2 * k1) + 2 * (u + dt / 2 * k2) + (u + dt * k3))
        u += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ts.append((k + 1) * dt); us.append(u); xs.append(x)
    return np.array(ts), np.array(us), np.array(xs)


def read(code, dt):
    D = np.genfromtxt(f"results/{code}/A_dt{dt}.csv", delimiter=",", names=True)
    return D["t"], D["ux"], D["x"]


tr, ur, xr = reference()
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

# 1. velocity
tn, un, xn = read("neofoam", DT_SHOW)
to, uo, xo = read("openfoam", DT_SHOW)
ax[0].plot(tr, ur, "k-", lw=1, label="reference (fine RK4)")
ax[0].plot(tn, un, "-", color="tab:blue", lw=2, label="NeoFOAM")
ax[0].plot(to, uo, "o", color="tab:orange", ms=4, mfc="none", markevery=4, label="OpenFOAM")
ax[0].set(xlabel="t [s]", ylabel="u [m/s]", title=f"Velocity (dt = {DT_SHOW} s)", xlim=(0, 0.1))
ax[0].legend()

# 2. position
ax[1].plot(tr, (xr - x0) * 1e3, "k-", lw=1, label="reference (fine RK4)")
ax[1].plot(tn, (xn - x0) * 1e3, "-", color="tab:blue", lw=2, label="NeoFOAM")
ax[1].plot(to, (xo - x0) * 1e3, "o", color="tab:orange", ms=4, mfc="none", markevery=4, label="OpenFOAM")
ax[1].set(xlabel="t [s]", ylabel="x − x0 [mm]", title=f"Position (dt = {DT_SHOW} s)")
ax[1].legend()

# 3. NeoFOAM vs OpenFOAM directly: the speed of both codes at the same states, every dt.
#    OpenFOAM is one step ahead (double first step of an injected particle), so the
#    OpenFOAM value at step k is paired with the NeoFOAM value at step k + 1.
largest = 0.0
for dt, c in zip(DTS, ["tab:purple", "tab:blue", "tab:green", "tab:orange", "tab:red"]):
    tn_, un_, _ = read("neofoam", dt)
    to_, uo_, _ = read("openfoam", dt)
    kn = np.round(tn_ / float(dt)).astype(int)
    ko = np.round(to_ / float(dt)).astype(int) + 1
    common, i_n, i_o = np.intersect1d(kn, ko, return_indices=True)
    largest = max(largest, np.abs(un_[i_n] - uo_[i_o]).max())
    ax[2].plot(uo_[i_o], un_[i_n], "o", ms=3, color=c, label=f"dt = {dt} s")
ax[2].plot([0, u0], [0, u0], "k-", lw=1, label="NeoFOAM = OpenFOAM")
ax[2].set(xlabel="u OpenFOAM [m/s]", ylabel="u NeoFOAM [m/s]",
          title=f"Same speed in both codes (largest difference {largest:.0e} m/s)", aspect="equal")
ax[2].title.set_fontsize(10)
ax[2].legend(fontsize=8)

fig.suptitle("Case A: one 50 µm particle slowed by drag (sphereDrag + Euler) — "
             "NeoFOAM and OpenFOAM give the same result", fontsize=12)
fig.tight_layout()
fig.savefig("figures/ref_A_simple.png", dpi=150)
print(f"largest |u_NeoFOAM - u_OpenFOAM| = {largest:.1e} m/s; written figures/ref_A_simple.png")
