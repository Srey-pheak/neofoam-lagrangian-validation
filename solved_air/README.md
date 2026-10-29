# Benchmark 2 – One-way coupled solver (live carrier)

This benchmark covers `neoIcoKinematicParcelFoam`: the incompressible PISO solve
of `neoIcoFoam` plus the Lagrangian particles moving in that solved flow. The
particles feel the fluid but do not act on it (one-way coupling).

It is separate from [Benchmark 1](../frozen_air/README.md), which
freezes the carrier and verifies the particle physics (Q3–Q8, Q1, Q2, Q10) in
isolation. Benchmark 2 verifies the **hand-over of the live flow** and the
coupled executable; it does not re-verify the particle physics.

> **“Coupled” here means “fluid and particles in one time loop”.** It is not
> two-way coupling: no momentum source is fed back to the fluid. See the code
> notes below.

## What is tested

| Test | Case | Reference | Runs | Claim |
|---|---|---|---|---|
| [05_uniform_flow](05_uniform_flow/README.md) | 1 m box, 10³ cells, solved uniform flow `(2 0 0)` | Test 3 (`C_RK2.csv`, frozen uniform flow) | 1 | The coupled solver reproduces the frozen-carrier result within round-off; the solved `U` is handed over every step |
| [06_cavity](06_cavity/README.md) | 1 × 1 × 0.1 m lid-driven cavity, 20 × 20 × 1 cells, air starts at rest and is solved every step | Python replay of the same steps on the `U` written each step, plus two wrong-air controls | 1 (smoke) | The particles are moved by the air **of the current step**, in the **correct cell**; a stale (`U(t_{n-1})`) or frozen (`U(0)`) air is detected |

**Status:** test 5 set up and PASS (max `|Δx| = 7.4e-15 m`, max `|Δu| = 1.1e-13 m/s`;
commit `21e102853`). Test 6 set up and smoke-PASS (correct air `2.2e-16 m`, one-step-old
`4.0e-3 m`, frozen `8.6e-1 m`); its recorded run still has to be made on a clean
commit (see [06_cavity/README.md](06_cavity/README.md#result)).

## What is not covered yet

| Not tested | Why it matters | Planned test |
|---|---|---|
| One-way enforced | Show that particles do not change the flow | G: run the same case with and without particles; `p`/`U` must match to round-off |
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

# Benchmark 2, test 5
cd neofoam-lagrangian-validation/solved_air/05_uniform_flow
scripts/run_E.sh

# Benchmark 2, test 6 (keeps the time folders the replay reads)
cd ../06_cavity
scripts/run_F.sh
```

For Benchmark 1, see its [README](../frozen_air/README.md#5-how-to-run).

## Layout

```
solved_air/
├── README.md                 this file
├── RECORD.md                 one line per run
├── draw_setups.py            initial-state schematics for tests 5 and 6
├── 05_uniform_flow/         Test 5: case, scripts, results, recorded run
│   ├── 0.orig/ constant/ system/ Allrun Allclean
│   ├── scripts/run_E.sh
│   ├── scripts/compare_to_C.py
│   ├── figures/E_setup.png
│   ├── results/
│   └── runs/
└── 06_cavity/                 Test 6: live non-uniform air (lid-driven cavity)
    ├── 0.orig/ constant/ system/ Allrun Allclean
    ├── scripts/run_F.sh
    ├── scripts/replay_F.py
    ├── figures/F_setup.png
    ├── results/
    └── runs/
```

## Related

- [Benchmark 1 – frozen carrier](../frozen_air/README.md)
- [Test 5 – setup, metric, result](05_uniform_flow/README.md)
- [Test 6 – setup, metric, result](06_cavity/README.md)
- Historical construction notes: [`../temp_coupled_step.md`](../temp_coupled_step.md)
