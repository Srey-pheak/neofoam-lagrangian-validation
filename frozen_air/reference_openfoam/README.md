# Reference – NeoFOAM vs OpenFOAM v2406

**Claim:** for the same case, NeoFOAM `neoIcoUncoupledKinematicParcelFoam` computes the same particle
motion as OpenFOAM v2406 `icoUncoupledKinematicParcelFoam`.

**Cases.** The Benchmark 1 cases (1, 2, 3), changed only where OpenFOAM lacks a model: OpenFOAM has
`sphereDrag` (Putnam) and `Euler` of our models, but no `stokesDrag`, `SchillerNaumannDrag` or RK2.
So every reference case uses **sphereDrag + Euler**, and NeoFOAM runs the identical case.

| Case | From | dt [s] | endTime [s] | Mesh | Particles |
|---|---|---|---|---|---|
| `A_dt<dt>` | 01_integrators | 2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4 | 0.2 | 1 cell | 1 (50 µm, u0 = 0.1 m/s), no gravity |
| `B` | 02_drag_settling | 5e-5 | 8 | 1 cell | 5 sizes, gravity |
| `C` | 03_motion_3d_exit | 5e-4 | 0.3 | 10×10×10 cells | 3, wind (2 0 0), gravity |

**Known differences between the codes**
- OpenFOAM moves a particle injected at SOI 0 **twice in its first time step**: its output at t_k
  equals the NeoFOAM state at t_k + dt. `compare.py` compares both alignments (offset 0 and 1 step).
- OpenFOAM **tracks through faces** and splits a step where a particle crosses a face (`maxCo`);
  NeoFOAM takes one step per dt. On a one-cell mesh (1, 2) this never happens; in test 3 it does, so
  small differences are expected there.
- OpenFOAM **removes** a particle that leaves; NeoFOAM keeps it as inactive.

**Comparisons (`compare.py`)**
1. NeoFOAM vs OpenFOAM, every case: max |Δu| / max|u| and max |Δx| for the best alignment
2. Test 1: error vs a fine RK4 reference of the same equation and the order of accuracy, both codes
3. Test 2: terminal velocity at t_end vs the force balance, both codes
4. Test 3: exit time of the particle that leaves through x = 1 m

**How to run**
```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/frozen_air/reference_openfoam
scripts/run_reference.sh                          # both codes, all cases (test 2 takes the longest)
scripts/run_reference.sh openfoam                 # OpenFOAM only;  ... neofoam  NeoFOAM only
scripts/run_reference.sh both A                   # one case group: A, A_dt1e-3, B or C
python3 scripts/compare.py | tee results/summary.txt
python3 scripts/plot_A_simple.py; python3 scripts/plot_B_simple.py; python3 scripts/plot_C_simple.py   # simple figures
```

**Folders**
- `runs/openfoam/<case>/`, `runs/neofoam/<case>/`: the run folders (separate per code; can be deleted)
- `results/openfoam/`, `results/neofoam/`: one CSV per case, same columns `t,id,x,y,z,ux,uy,uz,d`
- `scripts/`: `run_reference.sh`, `read_openfoam.py` (OpenFOAM `lagrangian/` output → CSV), `compare.py`
- `figures/`: `ref_A_simple.png` (test 1: velocity, position, speed point by point), `ref_A_convergence.png` (test 1: error vs dt), `ref_B_simple.png` (test 2: falling speed over time, final speed per size), `ref_C_simple.png` (test 3: paths from the side and from above, velocity in x, y, z)
