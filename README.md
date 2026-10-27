# NeoFOAM Lagrangian particle solver: validation cases

This repository holds only the **inputs and scripts** of the tests: case files,
run scripts, comparison and plot scripts. It holds no run output and no results.
Every number claimed in the `RECORD.md` files can be reproduced by running the
scripts against the solver commit given there.

Two benchmark suites for the NeoFOAM Lagrangian particle solver. They test
different programs and make different claims; keep them separate.

| Suite | Solver | Carrier | Scope | Status |
|---|---|---|---|---|
| [benchmark_uncoupleV1](benchmark_uncoupleV1/README.md) | `neoIcoUncoupledKinematicParcelFoam` | **frozen**, read once from `0/U` | particle physics: integrators, drag, gravity, cell location, exit (A–D, checks, CPU/GPU, OpenFOAM comparison) | A–D PASS |
| [benchmark_coupledV1](benchmark_coupledV1/README.md) | `neoIcoKinematicParcelFoam` | **live**, PISO solves `U` every step | hand-over of the solved flow to the particles; coupled executable (E) | E PASS |

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
cd neofoam-lagrangian-validation/benchmark_uncoupleV1
A_integrators/scripts/run_A.sh
B_forcesDrag/scripts/run_B.sh
C_3dBoundaries/scripts/run_C.sh
D_cellularFlow/scripts/run_D.sh

# Benchmark 2 (live carrier); compares to the output of test C, so run C first
cd neofoam-lagrangian-validation/benchmark_coupledV1/E_coupledUniform
scripts/run_E.sh
```

For the exact pass criteria and recorded results, see each suite's `README.md`
and `RECORD.md`.

## Layout

```
neofoam-lagrangian-validation/
├── README.md                     this index
├── benchmark_uncoupleV1/         Benchmark 1 — frozen carrier
│   ├── A_integrators/
│   ├── B_forcesDrag/
│   ├── C_3dBoundaries/
│   ├── D_cellularFlow/
│   ├── checks/
│   ├── cpu_vs_gpu/
│   ├── reference_openfoam/
│   ├── README.md
│   └── RECORD.md
├── benchmark_coupledV1/          Benchmark 2 — live carrier, one-way
│   ├── E_coupledUniform/
│   ├── README.md
│   └── RECORD.md
└── .gitignore                    keeps run output (runs/, results/, time folders) out of git
```

## Naming

- `benchmark_uncoupleV1` keeps its original name to avoid rewriting the paths in
  the existing Benchmark 1 scripts and results. It is the **uncoupled**
  (frozen-carrier) suite.
- `benchmark_coupledV1` is the **coupled** (live-carrier, one-way) suite.
- “Coupled” in the folder name means “fluid and particles in one time loop”, not
  two-way coupling. For two-way coupling, a momentum source from the particles
  would have to be added to the fluid equation; that is not implemented.
