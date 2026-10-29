# Benchmark 1 – run record

One line per test run, newest first. Solver = NeoFOAM git commit of the code that was built
("+ uncommitted" = the lagrangian code had changes not yet committed: commit before the next runs,
so every result points to an exact version). Details of each run: `<test>/results/summary.txt`,
data `<test>/results/*.csv`, figures `<test>/figures/`.

> **Commit hashes re-labelled on 2026-10-26.** The commit messages were rewritten (a
> `Co-Authored-By` trailer was removed); the trees are byte-identical, so every result below
> still refers to the same code. Old hashes: `cbb1efd9b` → `ff63558b8`,
> `3a38fbeaa` → `065cc52ff`.

| Date | Test | Runs | Solver (NeoFOAM) | Executor | Result | Note |
|---|---|---|---|---|---|---|
| 2026-10-20 | cpu_vs_gpu (1, 2, 3) | 3 + 3 | ff63558b8 + uncommitted (locateFindCell/locateWalk made public for CUDA) | Serial vs GPU (RTX 5060, CUDA 13.1) | **PASS** | test 1 bitwise identical; test 2 max difference 4.4e-16, test 3 2.8e-17 (round-off). GPU slower (test 2: 21.5 s vs 1.7 s): few particles, findCell copies GPU↔CPU every step. CPU results identical to the earlier 1, 2, 3 runs |
| 2026-10-20 | reference_openfoam 3 (rerun) | 1 + 1 | ff63558b8 + uncommitted; OpenFOAM v2406 | Serial | **same as OpenFOAM** | sphereDrag + Euler. Velocities differ ≤ 2.8e-5 m/s (OpenFOAM plotted one step ahead); exiting particle identical. New figure ref_C_simple.png |
| 2026-10-20 | 3 motion 3D and exit | 2 | ff63558b8 + uncommitted | Serial | **PASS** | RK2: u error ≤ 1.4e-4, x ≤ 5 µm, exit in the exact step. Euler: u error ≤ 9.8e-3, x ≤ 2 mm, exit 1 step late (limit 2). P2 frozen after leaving in both |
| 2026-10-20 | 2 drag and settling | 6 | ff63558b8 + uncommitted | Serial | **PASS** | stokes, sphere, SN, each RK2 and Euler. Stokes error RK2 ≤ 2.6e-4, Euler ≤ 1.1e-2 (10 µm). Terminal velocities ≤ 3e-13 for all. t95 matches theory (10 µm: 0.0027 vs 0.0023 only because output is every 0.002 s) |
| 2026-10-20 | reference_openfoam (1, 2, 3) | 7 + 7 | ff63558b8 + uncommitted; OpenFOAM v2406 | Serial | **same as OpenFOAM** | sphereDrag + Euler. Test 1: identical to ~1e-16 at every dt (offset 1 = OpenFOAM one step ahead), both codes order 1.00. Test 2: terminal velocities identical to 1e-13; transient differs only by the 1-step shift (output every 40 steps), then round-off. Test 3: exiting particle identical (3e-16); particles crossing faces differ 1e-5 / 6e-7 (OpenFOAM splits steps at faces) |
| 2026-10-20 | 1 integrators | 16 | ff63558b8 + uncommitted | Serial | **PASS** | Euler p = 1.00, RK2 p = 2.01 (finest dt). Time step: Euler slows down at dt = 0.001, 0.02, 0.05 s; RK2 speeds up at 0.05 s (2.6 τp) |
