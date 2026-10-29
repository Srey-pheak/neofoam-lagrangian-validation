# NeoFOAM Lagrangian particle solver: validation cases

This repository holds only the **inputs and scripts** of the tests: case files,
run scripts, comparison and plot scripts. It holds no run output and no results.
Every number claimed in the `RECORD.md` files can be reproduced by running the
scripts against the solver commit given there.

Two benchmark suites for the NeoFOAM Lagrangian particle solver. They test
different programs and make different claims; keep them separate.

| Suite | Solver | Carrier | Scope | Status |
|---|---|---|---|---|
| [frozen_air](frozen_air/README.md) | `neoIcoUncoupledKinematicParcelFoam` | **frozen**, read once from `0/U` | particle physics: integrators, drag, gravity, cell location, exit (1–4, checks, CPU/GPU, OpenFOAM comparison) | 1–4 PASS |
| [solved_air](solved_air/README.md) | `neoIcoKinematicParcelFoam` | **live**, PISO solves `U` every step | hand-over of the solved flow to the particles; coupled executable (5, 6) | 5 PASS; 6 PASS |

**Benchmark 1** verifies the particle model in isolation (Q3–Q8, Q1, Q2, Q10).
**Benchmark 2** verifies that the same model works with a carrier that is solved
every time step. Neither suite tests **two-way coupling**: the particles do not
feed momentum back to the fluid in either solver.

## Requirements

| | Version used |
|---|---|
| NeoFOAM (solver) | fork [`Srey-pheak/neofoam`](https://github.com/Srey-pheak/neofoam), branch `feature/lagrangian-particles`, commit `21e102853` (the commit of each run is in `RECORD.md`) |
| OpenFOAM | v2406 (`blockMesh`, run functions, and the reference solver `icoUncoupledKinematicParcelFoam`) |
| Python | 3.14.4, numpy 2.3.5, matplotlib 3.10.7 |

Build NeoFOAM, then tell the scripts where its binaries are:

```bash
export NEOFOAM_BIN=<NeoFOAM>/build/production/bin
export NEOFOAM_BUILD=<NeoFOAM>/build        # only for cpu_vs_gpu (needs production/ and production-cuda/)
```

## Quick start

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
export NEOFOAM_BIN=<NeoFOAM>/build/production/bin

# Benchmark 1 (frozen carrier)
cd neofoam-lagrangian-validation/frozen_air
01_integrators/scripts/run_A.sh
02_drag_settling/scripts/run_B.sh
03_motion_3d_exit/scripts/run_C.sh
04_cellular_flow/scripts/run_D.sh

# Benchmark 2 (live carrier); compares to the output of test 3, so run test 3 first
cd neofoam-lagrangian-validation/solved_air/05_uniform_flow
scripts/run_E.sh

# Test 6 (live non-uniform air; keeps the time folders its replay reads)
cd ../06_cavity
scripts/run_F.sh
```

For the exact pass criteria and recorded results, see each suite's `README.md`
and `RECORD.md`.

## Layout

```
neofoam-lagrangian-validation/
├── README.md                     this index
├── frozen_air/         Benchmark 1 — frozen carrier
│   ├── 01_integrators/
│   ├── 02_drag_settling/
│   ├── 03_motion_3d_exit/
│   ├── 04_cellular_flow/
│   ├── checks/
│   ├── cpu_vs_gpu/
│   ├── reference_openfoam/
│   ├── README.md
│   └── RECORD.md
├── solved_air/          Benchmark 2 — live carrier, one-way
│   ├── 05_uniform_flow/
│   ├── 06_cavity/
│   ├── draw_setups.py
│   ├── README.md
│   └── RECORD.md
└── .gitignore                    keeps run output (runs/, results/, time folders) out of git
```

## Naming

- `frozen_air` keeps its original name to avoid rewriting the paths in
  the existing Benchmark 1 scripts and results. It is the **uncoupled**
  (frozen-carrier) suite.
- `solved_air` is the **coupled** (live-carrier, one-way) suite.
- “Coupled” in the folder name means “fluid and particles in one time loop”, not
  two-way coupling. For two-way coupling, a momentum source from the particles
  would have to be added to the fluid equation; that is not implemented.
