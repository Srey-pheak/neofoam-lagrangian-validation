# Reference – test 6, NeoFOAM vs OpenFOAM v2406

**Claim:** for the lid-driven cavity of [test 6](../06_cavity/README.md) (air at rest, solved every
step, one-way particles):
1. NeoFOAM `neoIcoKinematicParcelFoam` solves the **same air** as OpenFOAM `icoFoam`, every step;
2. in that air, the NeoFOAM particles move **exactly as OpenFOAM's particle code** moves them, until
   a particle crosses a cell face, where the two codes treat the step differently (see below).

This complements the Python replay of test 6: the replay proves that the particles use the air of
the current step, but it reads that air from NeoFOAM itself. Here the air is checked against an
independent solver, and the particle motion against OpenFOAM's Lagrangian library.

**Case.** Test 6, changed only where needed, identically in every run:

| Change | Why |
|---|---|
| `sphereDrag` + `Euler` (test 6: `stokesDrag` + `RK2`) | OpenFOAM has no `stokesDrag` and no RK2 (as the [Benchmark 1 reference](../../frozen_air/reference_openfoam/README.md)) |
| linear solver tolerance `1e-12` (test 6: `1e-6` for p, `1e-5` for U) | NeoN and OpenFOAM normalise the residual differently. At the test 6 tolerances the two airs differ by the solver convergence, max `1.7e-4 m/s`, not by the method; at `1e-12` the method is compared |
| OpenFOAM only: `pFinal`/`UFinal` solvers, `PIMPLE` dictionary, laminar stress schemes | required by `icoFoam`/`kinematicParcelFoam` v2406; NeoFOAM reads the plain entries |

**Runs** (`scripts/run_reference.sh`)

| Run | Code | What |
|---|---|---|
| `neofoam/F` | `neoIcoKinematicParcelFoam` | air + particles; `U` written every step |
| `openfoam/F_fluid` | `icoFoam` | the air, same equations as `neoIcoFoam` |
| `openfoam/F_chain` | `icoUncoupledKinematicParcelFoam`, one launch per step (600) | OpenFOAM's particles in the **NeoFOAM air**: before the step from `t_k`, the NeoFOAM `U(t_k)` is copied into `t_k/U` |
| `openfoam/F` | `kinematicParcelFoam` | information only: OpenFOAM's own one-way solver |

Why the step-by-step chain and not `kinematicParcelFoam`:
- `kinematicParcelFoam` solves the air with the `pimpleFoam` momentum equation
  (`turbulence->divDevReff(U)`), not `icoFoam`'s. On this 20 × 20 mesh its air differs from
  `icoFoam` by up to `3.8e-2 m/s` (a separate `pimpleFoam` run gives the same air to the last digit,
  so the particles are not the cause: the cloud is one-way).
- It moves the particles **before** solving the air of the step (NeoFOAM: after). So its particles
  use the air one step older.

Both would mix into a particle comparison. The chain feeds OpenFOAM's particle code exactly the air
NeoFOAM used, so only the particle code differs.

**Alignment.** OpenFOAM moves a just-injected particle twice in its first step (see the Benchmark 1
reference). The chain's first step runs in the start air `U(0) = 0`, in which a particle at rest
does not move, so this has no effect, and OpenFOAM at `t_k+1` holds the NeoFOAM state at `t_k`.

**Known difference: the face crossing.** OpenFOAM tracks a particle face to face: where it crosses a
face during a step, it splits the step and uses the air of each cell for its part. NeoFOAM finds the
cell at the start of the step and uses that cell's air for the whole step. Before a particle's first
face crossing the two must agree to round-off; at a crossing they differ, and with cell values of U
(the air jumps from cell to cell) the difference grows at later crossings.

**Metric / Pass (`scripts/compare.py`)**
1. air: max |U_NeoFOAM − U_icoFoam| over all cells and steps `< 1e-9 m/s`
2. particles: max |x_NeoFOAM − x_OpenFOAM| before each particle's first face crossing `< 1e-12 m`;
   the difference after it is reported, not judged
3. `kinematicParcelFoam`: reported only

## Result

2026-10-28, `build/production` (rebuilt 2026-10-28 from `21e102853` plus uncommitted comment-only
changes), Serial:

| | Result | Limit | |
|---|---|---|---|
| Air, 601 steps | max \|ΔU\| = `4.1e-11 m/s` | `< 1e-9` | PASS |
| Particles before the first face crossing (t = 0.155 … 0.415 s) | max \|Δx\| = `6.7e-16 m` | `< 1e-12` | PASS |
| Particles after it, to t = 3 s | max \|Δx\| = `3.1e-2 m` (P0), `2.9e-3 … 1.3e-2 m` (others) | reported | |
| `kinematicParcelFoam` | air vs `icoFoam` `3.8e-2 m/s`; particles vs NeoFOAM `2.4e-2 m` | reported | |

The first face crossing changes the position by `1e-5 … 1e-3 m` (one step with a part of the
air jump of the crossed face); later crossings add to it in steps (figure `ref_F_simple.png`,
top right). The paths stay on top of each other at plot scale.

**Conclusion.** The air is the same as OpenFOAM's to round-off. The particle model is the same as
OpenFOAM's to round-off. The remaining particle difference is the cell-crossing treatment of the
locate step (NeoFOAM: start cell for the whole step; OpenFOAM: face-to-face tracking). It is a
property of the method, not an error of either code, and it grows with the size of the air jump
between cells, i.e. it would shrink with a finer mesh or a smooth (`cellPoint`) interpolation.

## How to run

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/solved_air/reference_openfoam
NEOFOAM_BIN=<NeoFOAM>/build/production/bin scripts/run_reference.sh    # about 3 min (the chain: 600 launches)
python3 scripts/compare.py | tee results/summary.txt
python3 scripts/plot_F_fluid.py; python3 scripts/plot_F_simple.py
```

## ParaView

Every run folder has a `.foam` file; open it with ParaView's built-in OpenFOAM reader:

| File | Shows |
|---|---|
| `runs/neofoam/F/neofoam_F.foam` | NeoFOAM air (`U`, `p`) and particles |
| `runs/openfoam/F_fluid/openfoam_F_fluid.foam` | OpenFOAM `icoFoam` air |
| `runs/openfoam/F_chain/openfoam_F_chain.foam` | OpenFOAM particles in the NeoFOAM air |
| `runs/openfoam/F/openfoam_F.foam` | OpenFOAM `kinematicParcelFoam` air and particles |

NeoFOAM writes its particles only to a CSV; `../write_lagrangian.py` (called by `run_reference.sh`)
also writes them as `<time>/lagrangian/kinematicCloud/{positions,U,d,origId}`, the layout
OpenFOAM writes, so the reader lists the cloud `kinematicCloud` for both codes (tick it under
*Mesh Regions* → `lagrangian/kinematicCloud`). Remember OpenFOAM's chain is one step behind:
OpenFOAM at `t + dt` = NeoFOAM at `t`.

## Folders

- `runs/neofoam/F/`, `runs/openfoam/{F_fluid,F_chain,F}/`: run folders, every step kept (can be deleted)
- `results/neofoam/F.csv`, `results/openfoam/{F_chain,F}.csv`: particles, columns `t,id,x,y,z,ux,uy,uz,d`;
  `results/summary.txt`: the output of `compare.py`
- `scripts/`: `run_reference.sh`, `read_openfoam.py` (OpenFOAM `lagrangian/` output → CSV, and cell
  values of `U`), `compare.py`, `plot_F_fluid.py`, `plot_F_simple.py`
- `figures/`: `ref_F_fluid.png` (air: centre-line profiles at t = 0.1, 0.5, 3 s and the largest
  difference every step), `ref_F_simple.png` (particles: paths, position difference against time
  with each particle's first face crossing, velocity in x and y)
