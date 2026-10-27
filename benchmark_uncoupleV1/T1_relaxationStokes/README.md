# T1 – Drag relaxation (Stokes drag, no gravity)


**Claim:** the drag force and the time integrators (Euler, RK2) are correct.


**Setup**
- Domain: box 0.1 × 0.1 × 0.1 m, 1 cell, all faces patch `sides` (type `patch`)
- Air at rest: U = (0 0 0), rhoInf = 1.2 kg/m³, nu = 1.5e-5 m²/s (mu = 1.8e-5 Pa·s)
- Particle: solid sphere, d = 50 µm, rho0 = 2500 kg/m³, at (0.05 0.05 0.05) m, U0 = (0.1 0 0) m/s
- Forces: `stokesDrag` only (CdRe = 24), no gravity
- endTime 0.2 s (≈ 10 τp)

One single cell: the particle never crosses a face (it moves only u0·τp = 1.93 mm).

**Reference (exact)**
- τp = rho0 d² / (18 mu) = 0.01929 s
- u(t) = u0 · exp(−t/τp)
- x(t) = x0 + u0 · τp · (1 − exp(−t/τp))

**Run:** Euler and RK2, dt = 1e-2, 5e-3, 2.5e-3, 1.25e-3, 6.25e-4, 3.125e-4.
NeoFOAM only: OpenFOAM has no `stokesDrag`.

**Metric:** e(dt) = |u_num − u_exact| / u0 at t = 0.2 s (and the same for x);
order p = log2(e(dt) / e(dt/2)).

**Pass:** Euler p ≈ 1, RK2 p ≈ 2 (±0.1).
