#!/usr/bin/env python3
"""CPU vs GPU: the same NeoFOAM case on the CPU (Serial) and on the GPU.

Reads results/<case>_CPU.csv and results/<case>_GPU.csv (t,id,x,y,z,ux,uy,uz,d) for A, B, C.
Prints per case: same number of rows, bitwise identical or not, largest difference, and the run
times (results/runtime.txt). PASS if the largest difference is < 1e-12 (round-off).

Writes figures/cpu_gpu.png: CPU = line, GPU = circles. Run from the cpu_vs_gpu folder:
    python3 scripts/plot_cpu_gpu.py
"""
import os

import numpy as np
import matplotlib.pyplot as plt

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
LIMIT = 1e-12

runtime = {}
if os.path.exists("results/runtime.txt"):
    for line in open("results/runtime.txt"):
        k, v = line.split()
        runtime[k] = float(v)


def read(name):
    f = f"results/{name}.csv"
    return np.genfromtxt(f, delimiter=",", names=True) if os.path.exists(f) else None


def by_particle(D):
    return {int(i): D[D["id"] == i] for i in np.unique(D["id"])}


os.makedirs("figures", exist_ok=True)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
ok_all = True
print(f"{'case':5s} {'rows CPU':>9s} {'rows GPU':>9s} {'bitwise':>8s} {'max |GPU - CPU|':>16s} "
      f"{'time CPU [s]':>13s} {'time GPU [s]':>13s}  verdict")
for a, case in zip(ax, "ABC"):
    C, G = read(f"{case}_CPU"), read(f"{case}_GPU")
    if C is None or G is None:
        print(f"{case:5s} missing: run scripts/run_cpu_gpu.sh {case}")
        ok_all = False
        continue
    same_shape = C.shape == G.shape
    bitwise = open(f"results/{case}_CPU.csv").read() == open(f"results/{case}_GPU.csv").read()
    diff = max(np.abs(C[k] - G[k]).max() for k in C.dtype.names) if same_shape else np.inf
    ok = same_shape and diff < LIMIT
    ok_all &= ok
    print(f"{case:5s} {len(C):9d} {len(G):9d} {str(bitwise):>8s} {diff:16.2e} "
          f"{runtime.get(case + '_CPU', np.nan):13.1f} {runtime.get(case + '_GPU', np.nan):13.1f}  "
          f"{'PASS' if ok else 'FAIL'}")

    Pc, Pg = by_particle(C), by_particle(G)
    for pid, c in zip(sorted(Pc), COLORS):
        p, q = Pc[pid], Pg.get(pid)
        lbl = f"d = {p['d'][0] * 1e6:.0f} µm"
        if case == "A":
            a.plot(p["t"], p["ux"], "-", color=c, lw=2, label="CPU")
            a.plot(q["t"], q["ux"], "o", color=c, ms=5, mfc="none", markevery=8, label="GPU")
        elif case == "B":
            a.plot(p["t"], -p["uz"], "-", color=c, lw=2, label=lbl)
            a.plot(q["t"], -q["uz"], "o", color=c, ms=4, mfc="none", markevery=25)
        else:
            a.plot(p["x"], p["z"], "-", color=c, lw=2, label=f"P{pid}, {lbl}")
            a.plot(q["x"], q["z"], "o", color=c, ms=4, mfc="none", markevery=20)
    word = "bitwise identical" if bitwise else f"max difference {diff:.1e}"
    if case == "A":
        a.set(xlabel="t [s]", ylabel="u_x [m/s]", title=f"A: particle slowing down — {word}")
    elif case == "B":
        a.set(xlabel="t [s]", ylabel="falling speed [m/s]", xlim=(0, 2.5), title=f"B: falling from rest — {word}")
        a.plot([], [], "k-", lw=2, label="CPU (lines)")
        a.plot([], [], "ko", ms=4, mfc="none", label="GPU (circles)")
    else:
        a.axvline(1.0, color="k", lw=1)
        a.set(xlabel="x [m]", ylabel="z [m]", title=f"C: side view in the wind — {word}")
        a.plot([], [], "k-", lw=2, label="CPU (lines)")
        a.plot([], [], "ko", ms=4, mfc="none", label="GPU (circles)")
    a.legend(fontsize=8)

fig.suptitle("NeoFOAM on the CPU (Serial) and on the GPU (RTX 5060): same cases, same result", fontsize=12)
fig.tight_layout()
fig.savefig("figures/cpu_gpu.png", dpi=150)
print(f"\nOVERALL: {'PASS' if ok_all else 'FAIL'}   (limit {LIMIT:g}; figures/cpu_gpu.png)")
