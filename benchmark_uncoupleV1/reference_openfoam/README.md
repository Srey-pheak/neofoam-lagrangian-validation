# Reference – NeoFOAM vs OpenFOAM v2406

**Claim:** for the same case, NeoFOAM `neoIcoUncoupledKinematicParcelFoam` computes the same particle
motion as OpenFOAM v2406 `icoUncoupledKinematicParcelFoam`.

**Cases.** The Benchmark 1 cases (A, B, C), changed only where OpenFOAM lacks a model: OpenFOAM has
`sphereDrag` (Putnam) and `Euler` of our models, but no `stokesDrag`, `SchillerNaumannDrag` or RK2.
So every reference case uses **sphereDrag + Euler**, and NeoFOAM runs the identical case.

| Case | From | dt [s] | endTime [s] | Mesh | Particles |
|---|---|---|---|---|---|
| `A_dt<dt>` | A_integrators | 2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4 | 0.2 | 1 cell | 1 (50 µm, u0 = 0.1 m/s), no gravity |
| `B` | B_forcesDrag | 5e-5 | 8 | 1 cell | 5 sizes, gravity |
| `C` | C_3dBoundaries | 5e-4 | 0.3 | 10×10×10 cells | 3, wind (2 0 0), gravity |

**Known differences between the codes**
- OpenFOAM moves a particle injected at SOI 0 **twice in its first time step**: its output at t_k
  equals the NeoFOAM state at t_k + dt. `compare.py` compares both alignments (offset 0 and 1 step).
- OpenFOAM **tracks through faces** and splits a step where a particle crosses a face (`maxCo`);
  NeoFOAM takes one step per dt. On a one-cell mesh (A, B) this never happens; in C it does, so
  small differences are expected there.
- OpenFOAM **removes** a particle that leaves; NeoFOAM keeps it as inactive.

**Comparisons (`compare.py`)**
1. NeoFOAM vs OpenFOAM, every case: max |Δu| / max|u| and max |Δx| for the best alignment
2. A: error vs a fine RK4 reference of the same equation and the order of accuracy, both codes
3. B: terminal velocity at t_end vs the force balance, both codes
4. C: exit time of the particle that leaves through x = 1 m

**How to run**
```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/benchmark_uncoupleV1/reference_openfoam
scripts/run_reference.sh                          # both codes, all cases (B takes the longest)
scripts/run_reference.sh openfoam                 # OpenFOAM only;  ... neofoam  NeoFOAM only
scripts/run_reference.sh both A                   # one case group: A, A_dt1e-3, B or C
python3 scripts/compare.py | tee results/summary.txt
python3 scripts/plot_A_simple.py; python3 scripts/plot_B_simple.py; python3 scripts/plot_C_simple.py   # simple figures
```

**Folders**
- `runs/openfoam/<case>/`, `runs/neofoam/<case>/`: the run folders (separate per code; can be deleted)
- `results/openfoam/`, `results/neofoam/`: one CSV per case, same columns `t,id,x,y,z,ux,uy,uz,d`
- `scripts/`: `run_reference.sh`, `read_openfoam.py` (OpenFOAM `lagrangian/` output → CSV), `compare.py`
- `figures/`: `ref_A_simple.png` (case A: velocity, position, speed point by point), `ref_A_convergence.png` (case A: error vs dt), `ref_B_simple.png` (case B: falling speed over time, final speed per size), `ref_C_simple.png` (case C: paths from the side and from above, velocity in x, y, z)
