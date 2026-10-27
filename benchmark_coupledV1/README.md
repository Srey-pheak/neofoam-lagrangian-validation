# Benchmark 2 – One-way coupled solver (live carrier)

This benchmark covers `neoIcoKinematicParcelFoam`: the incompressible PISO solve
of `neoIcoFoam` plus the Lagrangian particles moving in that solved flow. The
particles feel the fluid but do not act on it (one-way coupling).

It is separate from [Benchmark 1](../benchmark_uncoupleV1/README.md), which
freezes the carrier and verifies the particle physics (Q3–Q8, Q1, Q2, Q10) in
isolation. Benchmark 2 verifies the **hand-over of the live flow** and the
coupled executable; it does not re-verify the particle physics.

> **“Coupled” here means “fluid and particles in one time loop”.** It is not
> two-way coupling: no momentum source is fed back to the fluid. See the code
> notes below.

## What is tested

| Test | Case | Reference | Runs | Claim |
|---|---|---|---|---|
| [E_coupledUniform](E_coupledUniform/README.md) | 1 m box, 10³ cells, solved uniform flow `(2 0 0)` | Test C (`C_RK2.csv`, frozen uniform flow) | 1 | The coupled solver reproduces the frozen-carrier result within round-off; the solved `U` is handed over every step |

**Status:** E set up and PASS (max `|Δx| = 7.4e-15 m`, max `|Δu| = 1.1e-13 m/s`;
commit `21e102853`).

## What is not covered yet

| Not tested | Why it matters | Planned test |
|---|---|---|
| Live **non-uniform** `U` | E is uniform, so a stale or wrong-cell `U` would still pass. This is the main gap. | F: steady Couette or Poiseuille flow with an analytical `U(y)`, particle trajectories against an ODE reference |
| One-way enforced | Show that particles do not change the flow | G: run the same case with and without particles; `p`/`U` must match to round-off |
| Time-varying `U` | E reaches a steady state immediately | F/G with a ramp or a decaying flow |
| `locate walk` in the coupled solver | `CarrierField` is built once from the start `U`; walk would not see the flow change. The coupled program currently rejects walk. | rebuild the carrier field each step, or interpolate from `Ucells` in `locateWalk` |
| GPU | `locateFindCell` copies `Ucells` to the host every step | after a walk-based coupled locate |
| Two-way coupling | Not implemented; `coupled` in `kinematicCloudProperties` is ignored by NeoFOAM | separate solver work |

## How the code differs from Benchmark 1

The particle model is shared. The coupled program differs in **what it passes
to `ParticleSolver::locateAndInterpolate`**:

- Uncoupled: `Ufrozen`, copied once from `0/U` at `t = 0`
  ([`createFields.H`](../../NeoFOAM/examples/neoIcoUncoupledKinematicParcelFoam/createFields.H)).
- Coupled: `U.internalVector()` after the PISO loop
  ([`neoIcoKinematicParcelFoam.cpp`](../../NeoFOAM/examples/neoIcoKinematicParcelFoam/neoIcoKinematicParcelFoam.cpp)),
  i.e. the cell values of the flow just solved.

`ParticleSolver::locateAndInterpolate` takes that vector as its third argument
and `locateFindCell` uses `UV[celli]`. `ParticleSolver::evolve` is identical in
both programs and never writes back to the fluid.

## How to run

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc

# Benchmark 2, test E
cd neofoam-lagrangian-validation/benchmark_coupledV1/E_coupledUniform
scripts/run_E.sh
```

For Benchmark 1, see its [README](../benchmark_uncoupleV1/README.md#5-how-to-run).

## Layout

```
benchmark_coupledV1/
├── README.md                 this file
├── RECORD.md                 one line per run
└── E_coupledUniform/         Test E: case, scripts, results, recorded run
    ├── 0.orig/ constant/ system/ Allrun Allclean
    ├── scripts/run_E.sh
    ├── scripts/compare_to_C.py
    ├── results/
    ├── figures/
    └── runs/
```

## Related

- [Benchmark 1 – frozen carrier](../benchmark_uncoupleV1/README.md)
- [Test E – setup, metric, result](E_coupledUniform/README.md)
- Historical construction notes: [`../temp_coupled_step.md`](../temp_coupled_step.md)
