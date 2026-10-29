# Benchmark 1 – Verification of the particle physics (uniform carrier flow)

Benchmark 1 verifies the **current** `neoIcoUncoupledKinematicParcelFoam`: solid, spherical
particles (rho0 = 2500 kg/m³) in a **frozen, uniform** carrier flow, one-way coupling.
It proves that the solver computes the particle motion correctly (Q3–Q8), NeoFOAM only,
and that the basic mesh handling it has now (Q1, Q10) works. It does **not** test interpolation,
face tracking or periodic boundaries: the solver does not have them yet (see the end).

> **Scope (2026-10-26):** this folder is Benchmark 1 only — the frozen-carrier
> solver `neoIcoUncoupledKinematicParcelFoam`. The coupled case `05_uniform_flow`
> has moved to [`../solved_air/`](../solved_air/README.md)
> (Benchmark 2, live carrier). Tests 1–4 and the results here are unchanged.

## 1. What Benchmark 1 proves

| # | Claim | Proven by | Evidence (pass criterion) |
|---|---|---|---|
| 1 | Euler is 1st order in time | 1 | order p = 1 ± 0.1 between the two finest dt |
| 2 | RK2 is 2nd order in time | 1 | order p = 2 ± 0.1 between the two finest dt |
| 2b | A big time step: Euler still slows the particle down, RK2 speeds it up when dt > 2 τp | 1, part 2 | u(t) at dt = 0.001, 0.02, 0.05 s |
| 3 | The drag force is correct (size and direction) | 1, 2, 3 | match the exact Stokes solution (error < 1e-3 at fine dt) |
| 4 | Gravity with buoyancy is correct | 2, 3 | settling matches the exact Stokes solution for 5 sizes |
| 5 | Putnam drag (`sphereDrag`) is correct, Re 0.005 – 123 | 2 | terminal velocity = force balance, error < 1e-4 |
| 6 | Schiller–Naumann drag is correct, Re 0.005 – 123 | 2 | terminal velocity = force balance, error < 1e-4 |
| 7 | The motion is correct in all 3 components with moving air | 3 | exact per-component solution, error < 1e-3 (RK2), < 2e-2 (Euler) |
| 8 | Particles crossing cells are located without errors (Q1) | 3 | 10×10×10 mesh, trajectories match the exact solution (uniform air: a wrong cell would not show) |
| 8b | The correct cell is found in every step and its velocity is used (Q1, Q2), non-uniform air | 4 | bitwise equal to a reference that finds the cell by arithmetic |
| 8c | Particles are thrown out of the vortices (Stokes number effect) | 4 | positions at 1 and 2 τf; same as the smooth-flow reference |
| 9 | A particle leaving the domain is deactivated and stops (Q10) | 3 | exit within 2 dt of the exact crossing, frozen afterwards |
| 10 | Several particles in one cloud are independent (SoA indexing) | 2, 3 | each particle matches its own reference |
| 11 | Same result as OpenFOAM v2406 (sphereDrag + Euler) | reference_openfoam | identical to ~1e-16 for case 1; terminal velocities of test 2 to 1e-13 |
| 12 | The GPU kernel gives the same result on every executor | checks | Serial vs CPU bitwise identical; GPU < 1e-12 |
| 13 | Unsupported input stops the run with a clear error | checks | 6 / 6 broken cases stop with FATAL ERROR |

## 2. The tests

| Test | Case | Reference | Runs | Plots |
|---|---|---|---|---|
| [01_integrators](01_integrators/README.md) | box 0.1 m, 1 cell, still air, 1 particle (50 µm), Stokes drag | exact | 16 (Euler and RK2: 5 dt for the order, 3 dt for the size of the time step) | velocity, convergence, time step |
| [02_drag_settling](02_drag_settling/README.md) | column 0.5×0.5×35 m, 1 cell, still air, 5 particles (10 – 500 µm), drag + gravity | exact (Stokes), force balance | 6 (Stokes, Putnam, Schiller–Naumann; each RK2 and Euler) | falling speed over time, time to final speed, final speed vs size |
| [03_motion_3d_exit](03_motion_3d_exit/README.md) | box 1 m, 10×10×10 cells, wind (2 0 0), 3 particles, Stokes + gravity | exact per component | 2 (RK2 and Euler) | trajectories (x–z, x–y), error |
| [checks](checks/README.md) | case 3 (executors), case 1 broken 6 ways (input) | Serial run; FATAL ERROR | 2 – 3 + 6 | none (PASS/FAIL tables) |
| [04_cellular_flow](04_cellular_flow/README.md) | 2L × 2L, 64×64×1 cells, cellular flow (4 vortices), 4 × 100 particles (St 0.1 – 2), Stokes, no gravity | same-cell reference (bitwise), smooth-flow RK4 | 1 | positions at t = 0, 0.1, 0.2 s per St |
| [cpu_vs_gpu](cpu_vs_gpu/README.md) | reference cases 1, 2, 3 unchanged | CPU (Serial) run | 3 + 3 | CPU vs GPU: velocity, falling speed, side view |
| [reference_openfoam](reference_openfoam/README.md) | cases 1, 2, 3 with sphereDrag + Euler, in OpenFOAM v2406 and in NeoFOAM | OpenFOAM `icoUncoupledKinematicParcelFoam` | 7 per code | velocity, position, convergence of both codes (test 1) |

**Why the order of accuracy is tested only in test 1:** the integrator only sees `accel(u, x)`, never
the forces, so its order does not depend on the physics. Test 2 also runs every drag law with
Euler, to show both integrators reach the same terminal velocity.

> `T1_relaxationStokes/` is an earlier single-particle relaxation specification
> (no recorded results). The same physics is covered by test 1 (integrators) and test 2
> (drag); it is not part of the 1–4 status table.

## 3. What each part of the solver is tested by

| Solver part | File | 1 | 2 | 3 | checks |
|---|---|---|---|---|---|
| `StokesDrag` | `drag/StokesDrag.hpp` | ● | ● | ● | |
| `SphereDrag` (Putnam) | `drag/SphereDrag.hpp` | | ● | | |
| `SchillerNaumann` | `drag/SchillerNaumann.hpp` | | ● | | |
| `DragForce` (+ `Sp`) | `force/DragForce.hpp` | ● | ● | ● | |
| `GravityForce` | `force/GravityForce.hpp` | | ● | ● | |
| `EulerIntegrator` | `integrator/EulerIntegrator.hpp` | ● | ● | ● | ● |
| `RK2Integrator` | `integrator/RK2Integrator.hpp` | ● | ● | ● | ● |
| GPU kernel `ParticleSolver::evolve` | `src/lagrangian/ParticleSolver.cpp` | ● | ● | ● | ● |
| Q1/Q2 `locateFindCell` (findCell, cell value) | `src/lagrangian/ParticleSolver.cpp` | ○ | ○ | ● | ● |
| Q10 domain exit | `src/lagrangian/ParticleSolver.cpp` | | | ● | |
| Case reading and input checks | `createFields.H` | ○ | ○ | ○ | ● |

● tested directly · ○ runs but is not what the test checks

## 4. Not covered by Benchmark 1

| Not tested | Why | Planned |
|---|---|---|
| Interpolation (`cellPoint`) | needs `locate walk` (not verified); test 4 uses the cell value | test 4 with `locate walk` + `cellPoint` |
| Periodic (cyclic) boundaries | not with `findCell`: a leaving particle is deactivated (in test 4 most particles leave after a few τf) | test 4 with `locate walk` |
| Face-to-face tracking (Q9), GPU locate | `findCell` is used instead | MeshConnectivity on `rt.nfMesh`, checked against `findCell` |
| Walls (rebound, stick) | not implemented | if the thesis needs it |
| `analytical` integration | not implemented (input is rejected) | — |
| Dispersion, two-way coupling | not implemented | own tests when added; Benchmark 1 rerun as regression |
| Physical validity of the drag laws | validation (literature), not verification | — |

## 5. How to run

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/frozen_air
01_integrators/scripts/run_A.sh               && (cd 01_integrators    && python3 scripts/plot_A.py)
02_drag_settling/scripts/run_B.sh                && (cd 02_drag_settling     && python3 scripts/plot_B.py)
03_motion_3d_exit/scripts/run_C.sh              && (cd 03_motion_3d_exit   && python3 scripts/plot_C.py)
04_cellular_flow/scripts/run_D.sh              && (cd 04_cellular_flow   && python3 scripts/plot_D.py)
checks/scripts/check_executors.sh            # add GPU as argument on a GPU build
checks/scripts/check_input.sh
```
Every plot script prints a table and PASS or FAIL. Run the tests in order 1 → 2 → 3 → 4: each builds on
what the one before proved, so a failure points to one part of the solver.

**Layout of every test folder**
```
<test>/
├── README.md          claim, setup, reference, metric, pass
├── 0.orig/ constant/ system/ Allrun Allclean     reference case (never run here)
├── scripts/           run_*.sh (copies the case to runs/, runs, collects CSV), plot_*.py
├── runs/              one folder per run (can be deleted)
├── results/           CSV of every run (keep)
└── figures/           output of plot_*.py
```
All CSVs have the columns `t,id,x,y,z,ux,uy,uz,d` (id = injection order).

## 6. Status

| Test | Set up | Run | Result |
|---|---|---|---|
| 1 | ✅ | 16 | PASS: Euler p = 1.00, RK2 p = 2.01 |
| 2 | ✅ | 6 | PASS: terminal velocities ≤ 3e-13 (RK2 and Euler) |
| 3 | ✅ | 2 | PASS: RK2 error ≤ 1.4e-4, Euler ≤ 9.8e-3; exit within 1 step |
| 4 | ✅ | — (test run in scratchpad: PASS, difference 0.0) | — |
| cpu_vs_gpu | ✅ | 3 + 3 | PASS: 1 bitwise identical, 2 4.4e-16, 3 2.8e-17 |
| checks | ✅ | — | — |

Plot scripts were checked on synthetic exact-solution data (all PASS) before any real run.
