# Test E — Coupled uniform flow (consistency with Benchmark 1, Test C)

**Claim:** with a uniform carrier, the one-way coupled solver
`neoIcoKinematicParcelFoam` (PISO + particles) reproduces the frozen-carrier
trajectories of [Test C](../../benchmark_uncoupleV1/C_3dBoundaries/README.md)
within solver round-off. This verifies that the coupled program runs, that the
solved `U` is handed to the particle solver, and that the flow solve does not
perturb a uniform field. Any error in handing the air over changes the paths.

**What it does not prove:** that the particles see the *current* `U`. The
carrier is uniform and steady, so the result would be the same if `U` were
frozen at `t = 0`, and every cell has the same value. A non-uniform or
time-changing live-flow test is still needed (see Benchmark 2 README).

## Setup

- Domain: 1 m cube, 10 × 10 × 10 cells (`system/blockMeshDict`).
- Boundaries: inlet `U = (2 0 0)`, outlet `p = 0`, `sides` symmetry. With these
  boundaries, uniform flow is an exact solution of the incompressible equations.
- Fluid: `nu` and `rhoInf` from `constant/transportProperties`; PISO from
  `system/fvSchemes` and `system/fvSolution`.
- Particles: 3 solid spheres, `rho0 = 2500 kg/m³`:
  | id | d [µm] | initial U [m/s] | note |
  |---|---|---|---|
  | 0 | 50 | (0 0 0) | carried by the wind |
  | 1 | 40 | (0 1 0) | pushed sideways, then picked up |
  | 2 | 60 | (0 0 0) | leaves the domain |
- Forces: `stokesDrag` + `gravity`; integration RK2; `dt = 5e-4 s`;
  `endTime = 0.3 s`; output every step.
- `solution/locate findCell`; `interpolationSchemes/U cell`.

The boundary split and the case construction are described step by step in
[`../../temp_coupled_step.md`](../../temp_coupled_step.md) (historical notes;
this folder is the benchmark form).

## Reference

[`C_RK2.csv`](../../benchmark_uncoupleV1/C_3dBoundaries/results/C_RK2.csv):
the same three particles, same box and same `dt`, but with the carrier frozen at
`U = (2 0 0)` in `neoIcoUncoupledKinematicParcelFoam`.

## Metric / pass

- Same number of rows, same `t` column, same `id` column.
- `max |Δx|`, `max |Δy|`, `max |Δz| < 1e-9 m`.
- `max |Δux|`, `max |Δuy|`, `max |Δuz| < 1e-9 m/s`.
- `Δd = 0` exactly.

The limits are far below any physical effect; the observed difference is
round-off. The comparison is in `scripts/compare_to_C.py`.

## Result

Run `E_uniform`, 2026-10-26 11:29, solver commit `21e102853`, Serial:

| quantity | max \|E − C\| | limit | result |
|---|---|---|---|---|
| x | 7.4e-15 m | 1e-9 m | PASS |
| ux | 1.1e-13 m/s | 1e-9 m/s | PASS |
| y, z, uy, uz, d | 0 | — | PASS |

`y` and `z` are bitwise identical because the solved `u_y` and `u_z` are exactly
zero; only `x` and `u_x` show round-off.

1803 rows, `t = 0 … 0.3 s`. Full output: [`results/E_uniform.summary.txt`](results/E_uniform.summary.txt),
data: [`results/E_uniform.csv`](results/E_uniform.csv), log:
[`runs/E_uniform/log.neoIcoKinematicParcelFoam`](runs/E_uniform/log.neoIcoKinematicParcelFoam).

## How to run

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/benchmark_coupledV1/E_coupledUniform
scripts/run_E.sh                 # default run name: E_uniform
# or: scripts/run_E.sh E_uniform_rerun
```

`run_E.sh` copies the clean case to `runs/<name>/`, runs `blockMesh` and the
solver, copies the particle CSV to `results/<name>.csv`, and runs the comparison.
Set `NEOFOAM_BIN` if the production binary is somewhere else. `Allrun` removes
the time folders at the end; only `0/`, the logs and the particle CSV remain.

## Layout

```
E_coupledUniform/
├── README.md
├── 0.orig/ constant/ system/     clean case (no mesh, no time folders)
├── Allrun Allclean
├── scripts/run_E.sh              case -> runs/<name>, run, collect CSV
├── scripts/compare_to_C.py       compare with C_RK2.csv, PASS/FAIL
├── results/E_uniform.csv         particle data of the recorded run
├── results/E_uniform.summary.txt  PASS/FAIL table of the recorded run
├── figures/                      (plots, when added)
└── runs/E_uniform/               recorded run: mesh, 0/, logs, postProcessing/
```

## Note on `coupled false`

`constant/kinematicCloudProperties` contains `coupled false;`. This is an
OpenFOAM cloud keyword; the NeoFOAM reader does not read it (it parses
`solution/*` and `subModels/*` only). The one-way behaviour is enforced by the
code — `ParticleSolver::evolve` never writes back to the fluid — not by this
flag.
