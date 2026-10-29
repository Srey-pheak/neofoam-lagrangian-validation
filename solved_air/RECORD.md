# Benchmark 2 – run record

One line per test run, newest first. Solver = NeoFOAM git commit of the code that
was built (“+ uncommitted” = the lagrangian code had changes not yet committed:
commit before the next runs, so every result points to an exact version).
Details in `<test>/results/<run>.summary.txt`, data `<test>/results/*.csv`, logs
`<test>/runs/<run>/log.*`.

| Date | Test | Runs | Solver (NeoFOAM) | Executor | Result | Note |
|---|---|---|---|---|---|---|
| 2026-10-26 | 05_uniform_flow | 1 | 21e102853 | Serial | **PASS** | max `\|Δx\| = 7.4e-15 m`, max `\|Δu\| = 1.1e-13 m/s` vs `C_RK2.csv`; 1803 rows, `t = 0 … 0.3 s`; 1 particle leaves, as in C. Run 11:29, log in `05_uniform_flow/runs/E_uniform/` |

**Notes**

- The test 5 run was made before the folder split (`05_uniform_flow` was still under
  `frozen_air/`); the log therefore prints the old case path. The data
  and the code tree (`21e102853`) are the same.
- 1/2/3/4 were not rerun; Benchmark 1 covers them. Test 5 is the only coupled test
  so far.
- The comparison script is `05_uniform_flow/scripts/compare_to_C.py`; rerun it any
  time with `python3 scripts/compare_to_C.py`.

## Set up, awaiting a recorded run

| Date | Test | What | Why it is not recorded yet |
|---|---|---|---|
| 2026-10-28 | [06_cavity](06_cavity/README.md) | lid-driven cavity, live non-uniform and time-changing `U`; Python replay + two wrong-air controls | case and script are ready and a smoke run passed (correct air `2.2e-16 m`, one-step-old `4.0e-3 m`, frozen `8.6e-1 m`), but the smoke run used a tree with uncommitted changes. Rebuild from a clean commit, run `scripts/run_F.sh`, then move the row into the table above. |

## Planned

| Date | Test | What | Why |
|---|---|---|---|
| — | G | same case with and without particles | prove one-way coupling: `p`/`U` unchanged by the particles |
