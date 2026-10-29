# Test 4 – Cellular flow (Maxey 1987): non-uniform air, correct cell, Stokes number

**Claim:**
1. the particle is located in the **correct cell** in every step (Q1) and gets the velocity of
   that cell (Q2): in a non-uniform flow a wrong cell would give a wrong velocity
   (tests 1–3 have uniform air, so they cannot show this)
2. the particles behave as the Stokes number predicts: they are thrown out of the vortices

**Setup** (reference case in this folder; mesh, flow formula, particles and time step come from
`scripts/params.py`, written into the case by `scripts/make_case.py`)

- Domain: 2L × 2L × dz, L = 0.01 m, **64 × 64 × 1 cells** (32 per vortex, 1 cell thick in z),
  patch `sides` (type `patch`)
- Air: frozen cellular flow, 2 × 2 = 4 vortices, U0 = 0.1 m/s, τf = L / U0 = 0.1 s:
  - u_x =  U0 sin(πx/L) cos(πy/L)
  - u_y = −U0 cos(πx/L) sin(πy/L)

  The formula is in `system/setExprFieldsDict`. `Allrun` runs `setExprFields`, which evaluates it at
  every cell centre and writes `0/U` (one value per cell). The solver uses the cell value (`locate findCell`).
- 400 solid particles (rho0 = 2500 kg/m³), released at rest on a 10 × 10 grid (shifted by 0.3 cell,
  so no particle starts on a cell face), 100 per size:

| St = τp / τf | d [µm] | τp [s] |
|---|---|---|
| 0.1 | 36.0 | 0.01 |
| 0.5 | 80.5 | 0.05 |
| 1 | 113.8 | 0.1 |
| 2 | 161.0 | 0.2 |

  The diameter sets St: τp = ρp d² / (18 μ). Largest Re ≈ 1, so Stokes drag is used.
- Forces: `stokesDrag` only (no gravity); RK2, dt = 1e-4 s, endTime 1 s (10 τf), output every 50 steps

**References** (`scripts/plot_D.py`, independent of NeoFOAM)
1. **same-cell reference:** the same steps as the solver (find the cell, take its value of U, one RK2
   step), but the cell comes from arithmetic, cell = floor(x/dx) + nx·floor(y/dx), not from `findCell`.
   It uses the cell values of the run (`results/U_cells`, the `0/U` written by `setExprFields`).
   If NeoFOAM finds the right cell in every step, both agree to round-off.
2. **smooth-flow reference:** fine RK4 with the exact (smooth) flow. Shows the effect of using the
   cell value of U. Information only.

**Metric / Pass**
- max |x_NeoFOAM − x_same-cell| over all particles and outputs < 1e-9 m (cell size 3.1e-4 m)
- the same particles leave the domain in NeoFOAM and in the same-cell reference

**Note on leaving:** the domain has no periodic boundaries. The outer vortex boundaries are the domain
edges, so particles thrown out of the vortices reach an edge and leave (deactivated). After 1 s almost
all have left; the figure therefore shows t = 0, 0.1 and 0.2 s.

**How to run**
```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/frozen_air/04_cellular_flow
python3 scripts/make_case.py            # only after changing scripts/params.py
scripts/run_D.sh                        # 1 run, about 30 s -> results/D.csv, results/U_cells
python3 scripts/plot_D.py | tee results/summary.txt   # tables + PASS/FAIL, figures/D_particles.png
```

**Figure:** `figures/D_particles.png`: one row per St, positions at t = 0, 0.1 and 0.2 s on the streamlines.

**Folders:** reference case (`0.orig/ constant/ system/ Allrun Allclean`), `runs/`, `results/`, `scripts/`, `figures/`
