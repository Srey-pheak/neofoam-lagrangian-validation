# Benchmark 2 – run record

One line per test run, newest first. Solver = NeoFOAM git commit of the code that
was built (“+ uncommitted” = the lagrangian code had changes not yet committed:
commit before the next runs, so every result points to an exact version).
Details in `<test>/results/<run>.summary.txt`, data `<test>/results/*.csv`, logs
`<test>/runs/<run>/log.*`.

| Date | Test | Runs | Solver (NeoFOAM) | Executor | Result | Note |
|---|---|---|---|---|---|---|
| 2026-10-26 | E_coupledUniform | 1 | 21e102853 | Serial | **PASS** | max `\|Δx\| = 7.4e-15 m`, max `\|Δu\| = 1.1e-13 m/s` vs `C_RK2.csv`; 1803 rows, `t = 0 … 0.3 s`; 1 particle leaves, as in C. Run 11:29, log in `E_coupledUniform/runs/E_uniform/` |

**Notes**

- The E run was made before the folder split (`E_coupledUniform` was still under
  `benchmark_uncoupleV1/`); the log therefore prints the old case path. The data
  and the code tree (`21e102853`) are the same.
- A/B/C/D were not rerun; Benchmark 1 covers them. E is the only coupled test
  so far.
- The comparison script is `E_coupledUniform/scripts/compare_to_C.py`; rerun it any
  time with `python3 scripts/compare_to_C.py`.

## Planned

| Date | Test | What | Why |
|---|---|---|---|
| — | F | steady Couette/Poiseuille, non-uniform live `U` | E cannot detect a wrong-cell or stale `U` because the carrier is uniform |
| — | G | same case with and without particles | prove one-way coupling: `p`/`U` unchanged by the particles |
