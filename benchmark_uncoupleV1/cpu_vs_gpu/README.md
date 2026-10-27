# CPU vs GPU – the same NeoFOAM case on both

**Claim:** the solver gives the same result on the GPU as on the CPU.

**Runs** (`scripts/run_cpu_gpu.sh`, NeoFOAM only): the reference cases of A, B and C, unchanged, each twice

| Case | From | What | Steps |
|---|---|---|---|
| `A` | A_integrators | 1 particle slowing down, Stokes drag, Euler, dt = 1e-3 | 200 |
| `B` | B_forcesDrag | 5 particles falling from rest, sphereDrag + RK2, dt = 5e-5 | 160 000 |
| `C` | C_3dBoundaries | 3 particles in the wind, 1000 cells, Stokes + RK2, dt = 5e-4 | 600 |

| Run | Build | Executor |
|---|---|---|
| `<case>_CPU` | `NeoFOAM/build/production` | `-executor Serial` (as all other tests) |
| `<case>_GPU` | `NeoFOAM/build/production-cuda` | `-executor GPU` (RTX 5060) |

**Metric / Pass:** largest difference of every CSV column, GPU vs CPU, < 1e-12 (round-off).
Also printed: bitwise identical or not, and the run time of each run.

**How to run**
```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/benchmark_uncoupleV1/cpu_vs_gpu
scripts/run_cpu_gpu.sh                    # all cases (or: scripts/run_cpu_gpu.sh C)
python3 scripts/plot_cpu_gpu.py | tee results/summary.txt   # table + PASS/FAIL, figures/cpu_gpu.png
```
Needs the CUDA build (`build/production-cuda`, built with `-j 4`).

**Figure:** `figures/cpu_gpu.png`: A velocity, B falling speed, C side view; CPU = line, GPU = circles.

**Result (2026-10-20): PASS**

| Case | GPU vs CPU | Time CPU | Time GPU |
|---|---|---|---|
| A | bitwise identical | 0.6 s | 1.0 s |
| B | max difference 4.4e-16 | 1.7 s | 21.5 s |
| C | max difference 2.8e-17 | 0.6 s | 1.6 s |

The tiny differences in B and C are round-off (most likely: the GPU does some multiply-adds in one operation
and rounds once instead of twice). The GPU is slower here because there are only 1–5 particles
(nothing to share among thousands of GPU threads) and the default `locate findCell` copies the
particles GPU → CPU → GPU every step (B: 160 000 steps).
