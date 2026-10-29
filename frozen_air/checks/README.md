# Checks – executors and input validation (no plots)

## 1. Executors: `scripts/check_executors.sh [GPU]`

**Claim:** the particle kernel gives the same result on every executor.

Runs the test 3 case (3 particles, 1000 cells, crossing cells, one leaving the domain) with
`-executor Serial`, `CPU` and, with the argument `GPU`, `GPU`, and compares the CSVs with Serial.

**Pass:** Serial vs CPU bitwise identical; GPU identical or round-off (< 1e-12).

## 2. Bad input: `scripts/check_input.sh`

**Claim:** case input the solver does not support stops the run with a clear error.

Starts from the test 1 case and breaks one thing per run:

| Run | Change | Expected |
|---|---|---|
| `two_drag_models` | `stokesDrag;` + `sphereDrag;` | FATAL ERROR: only one drag model |
| `no_drag_model` | drag line removed | FATAL ERROR: a drag model is required |
| `unknown_force` | + `pressureGradient;` | FATAL ERROR: not supported |
| `analytical_scheme` | `U analytical;` | FATAL ERROR: not implemented |
| `SOI_not_zero` | `SOI 0.1;` | FATAL ERROR: only SOI 0 |
| `injection_outside_mesh` | position (5 5 5) | FATAL ERROR: outside the mesh |

**Pass:** 6 / 6 rejected. The messages are written to `results/check_input.txt`.

**How to run**
```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd neofoam-lagrangian-validation/frozen_air/checks
scripts/check_executors.sh
scripts/check_input.sh
```
