# particle_tracking.py

3D full-orbit test particles in saved formation runs
(`plasma_bench/particle_tracking.py`). This is step 1 toward the 3D model
(frc-ibl.12).

```bash
.venv/bin/python -m examples.particle_tracking verify outputs/particle_tracking/verify
.venv/bin/python -m examples.particle_tracking ions outputs/fine_uniform_h2 outputs/particle_tracking/ions_fine_uniform_h2
.venv/bin/python -m examples.particle_tracking rmf outputs/fine_uniform_h2 outputs/particle_tracking/rmf_even
.venv/bin/python -m examples.particle_tracking multipoles outputs/fine_uniform_h2 outputs/particle_tracking/multipoles_fine_uniform_h2
```

Options are `key=value`:
- `ions`: `start_us=10 end_us=30 count=4000 seed=1 window_us=2 trails=300 omega_dt=0.1`
- `rmf`: `start_us=10 duration_us=2 count=300 frequency=1e6 amplitude=5e-3 parity=even radius=0.16 half_length=0.6 omega_dt=0.3`
- `multipoles`: `start_us=10 end_us=30 count=4000 seed=1 cases=none,q3,q5,q7,q9,o7 omega_dt=0.1 workers=0` (q = 4-bar quadrupole, o = 8-bar octupole, number = % of B_e at the wall; one process per case unless `workers` is set)

## Fields and approximations
**Flux (`Run`).** The run saves the total flux on the gas nodes. It is
extended exactly as in `formation_display.py`:
1. The coil flux of the saved section vacuum currents is subtracted.
2. `CouplingMesh.harmonic_extension` fills the vacuum: outside the wall,
   inside an inner tube, and beyond the ends.
3. The coil flux is added back.

The gas-node flux is reproduced to 4e-17. Between the 0.1 µs saved states the
flux is linear in time.

**B = ∇ψ×θ̂/r.** It comes from the spline ψ = r² g. Against the display's
`np.gradient` field on the interior gas nodes it differs by 0.26% rms and
2.9% at most, which is the size of the grid discretization.

**E_θ = −(∂ψ/∂t)/r.** This is Faraday's law for the interpolated flux.
- **No −v×B term is added.** The single-fluid model has no toroidal field or
  flow, so −v×B = (v_r B_z − v_z B_r) θ̂ is purely azimuthal and is already
  part of the Faraday E_θ. Adding it would double count. The model's
  E_r = E_z = 0.
- **Consistency check.** Inside the separatrix at 10.1 µs, the Faraday E_θ
  (a central difference about the saved state) correlates 0.97 with the
  saved −(v×B)_θ. The rms difference is 28%: the resistive (Spitzer +
  Chodura) and numerical parts of Ohm's law.

**Coils.** Coil fields are vacuum Biot–Savart superposed on the saved field.
- The quadrupole is static (E = 0).
- The RMF antenna's E is −(dI/dt)A.
- Neither the plasma nor the sheet coil screens them in this model.

**Loading (`Run.load`).** Particles are placed inside the separatrix of the
saved state at `start_us`.
- **Cells:** probability ∝ n r, uniform in volume within a cell.
- **Positions:** kept only where the spline ψ < ψ_separatrix and between the
  closed-region ends.
- **Velocities:** a local Maxwellian at T = p/(2n) (Te = Ti), shifted by the
  cell's saved (v_r, v_z).
- **Not an equilibrium:** a local Maxwellian is not a function of the orbit
  invariants, so the population phase-mixes over the first few orbits.

**Losses.** Particles are lost at r = 0.1524 m, at the inner tube
r = 0.0381 m (annular runs) and at |z| = 1.1 m.

**Time step.** ω_c,max dt = 0.1 for ions (dt = 4.0–4.2 ns), with an integer
number of steps per saved interval.

## Verification (`verify`, 2026-10-20, outputs/particle_tracking/verify)
| Check | Reference | Result |
|---|---|---|
| Gyro-frequency, uniform 0.2 T (100 eV H⁺, 100 eV e⁻, 500 keV e⁻ with γ = 1.98) | qB/(γm) | Error −3.313e-3 / −8.321e-4 / −2.083e-4 at ωdt = 0.2 / 0.1 / 0.05. This equals the Boris angle 2 atan(ωdt/2) to 1e-15: second order. |
| Gyro-radius (circle fit) | u/(γω) | Exact to ≤ 2e-13 at every dt. The DKD radius is exact. |
| Energy, uniform B | — | ≤ 9e-15 relative after 20 turns. |
| E×B drift, E = 10 kV/m, B = 0.2 T | E×B/B² = −5×10⁴ m/s | ≤ 7e-9 relative (H⁺ and e⁻, ωdt = 0.5 and 0.1). |
| ∇B drift, B = 0.2(1 + x/1 m) T, 52 eV H⁺ (v⊥ = 10⁵ m/s, ρ/L = 5.2e-3) | v⊥²/(2ω_c L)·(B0/B_gc)² = 258.29 m/s | 258.30 m/s at every ωdt from 0.4 to 0.05. The 5.5e-5 residual is the O((ρ/L)²) finite-Larmor term. |
| Mirror loss cone: two loops, exact elliptic-integral flux passed through `FluxField`; R = 3.02, ρ/L = 0.015 | Isotropic 1 − √(1 − 1/R) = 18.2%; per particle, lost ⇔ sin²α < B₀/B_max | 17.75% lost (2000 protons). **0 of 2000** differ from the adiabatic per-particle prediction. |
| Quadrupole at the midplane wall | four infinite line currents | Ratio 1.0003. Azimuthal ripple ×2.0 at r_w/a = 0.76, from the (r/a)⁴ harmonic. 3.14 mT/kA. |
| Dossier quadrupole current (7% of 0.35 T at the wall) | BUILD_DOSSIER §6: ~8 kA | 7.8 kA |
| RMF antenna, Faraday's law | ∇×E = −∂B/∂t by central differences | 9e-7 relative residual |
| Energy, frozen FRC (10 µs) + 3D quadrupole, 200 H⁺ for 4 µs | exact | ≤ 1.5e-14 relative, ωdt = 0.2–0.025 |
| p_θ, frozen axisymmetric FRC, 200 H⁺ for 4 µs | exact invariant | max \|Δp_θ\|/(eψ_O) = 1.2e-3, 3.3e-4, 8.6e-5, 2.2e-5: **order 1.98, 2.00, 2.00** |
| p_θ, time-dependent FRC (10–14 µs, ψ linear in t between saved states) | exact invariant | 9.3e-4, 2.6e-4, 6.9e-5, 1.75e-5: **order 1.92, 2.00, 2.00** |

## Ion orbits in the formed FRC (`ions`, outputs/particle_tracking/ions_*)
**Setup.** 4000 thermal ions (H⁺, from the H2 runs), loaded at 10 µs and
tracked to 30 µs, in three field cases:
- no quadrupole;
- the straight quadrupole with 7% of B_e at the wall;
- 7% of B_e at the separatrix (the dossier's budgeted level).

B_e is |B_z| at the midplane wall at 10 µs. The quadrupole level is the mean
|B_quad| around that radius. All cases share the same initial particles.

**Orbit classes.** Classified over the first 2 µs (≈ 8 gyro-periods at B_e):
- **Betatron:** crosses the field null (B_z changes sign at least twice)
  and θ̇ keeps one sign, so the orbit encircles the axis. These are the
  meandering figure-8 orbits.
- **Mixed:** crosses the null, but θ̇ reverses.
- **Drift:** does not cross the null, and the magnetic moment varies by
  more than ×1.5.
- **Cyclotron:** does not cross the null, with μ_max/μ_min ≤ 1.5.
- **Lost in window:** lost within the first 2 µs.

The thresholds hardly matter: with ×1.2 or ×2 instead of ×1.5, the
cyclotron fraction stays ≤ 2%.

| Run (at 10 µs) | Option A, `fine_uniform_h2` | Inner tube, `fine_inner_tube` | Flux core, `fine_flux_core` |
|---|---|---|---|
| r_s / l_s / φ / Te+Ti | 11.4 cm / 0.86 m / 2.0 mWb / 157 eV | 10.4 cm (inner 4.2 cm) / 0.64 m / 1.24 mWb / 172 eV | 10.7 cm (inner 3.9 cm) / 0.61 m / 1.42 mWb / 176 eV |
| B_e, ρ_i at B_e | 0.245 T, 5.2 mm | 0.233 T, 5.8 mm | 0.239 T, 5.7 mm |
| s̄ = r_s/(4ρ_i,e) (repository convention) | 5.4 | 4.5 | 4.7 |
| s̄ with the local \|B\|, null → separatrix | 2.1 | 1.1 | 1.3 |
| Betatron / mixed / drift / cyclotron / lost in 2 µs | 13 / 46 / 41 / 0.3 / 0% | 19 / 40 / 26 / 0 / 15% | 20 / 42 / 28 / 0.1 / 10% |
| Confined at 20 / 30 µs, no quadrupole | 96.3 / 95.0% | 78.3 / 76.6% | 87.2 / 87.1% |
| Same, quadrupole 7% B_e at the wall (bar current) | 90.9 / 71.4% (5.5 kA) | 60.2 / 37.6% (5.2 kA) | 71.7 / 51.3% (5.3 kA) |
| Same, quadrupole 7% B_e at the separatrix | 87.0 / 62.1% (7.5 kA) | 51.6 / 27.4% (7.8 kA) | 62.4 / 38.3% (7.8 kA) |
| Where losses go, no quadrupole | ends 3.8%, wall 1.2% | inner tube 23.3% (median 10.9 µs) | inner tube 12.8% (median 10.5 µs) |
| Where losses go, quadrupole at the wall | outer wall 25.4%, ends 3.3% | inner tube 48.8%, outer wall 13.6% | inner tube 30.7%, outer wall 18.0% |

**Readings.**

**Ions are kinetic in this FRC.**
- Within 2 µs, 58–62% of the ions cross the field null.
- Fewer than 0.3% keep their magnetic moment (cyclotron class). The median
  μ_max/μ_min over 2 µs is 50–60.
- The repository's s̄ (ρ_i taken at B_e) is 4.5–5.4. Measured with the local
  field it is 1.1–2.1, because |B| falls to zero at the null (R ≈ 8 cm),
  where the gyroradius diverges. Only ~2 local gyroradii separate the null
  from the separatrix.

**Axisymmetric confinement.**
- In the axisymmetric fields, betatron ions are all confined: p_θ
  conservation forbids their escape.
- The ions lost come from the drift and mixed classes near the separatrix.
- The few cyclotron-class ions (10 of 4000) are mostly lost: adiabatic ions
  follow field lines out once the shrinking separatrix leaves them outside.

**Quadrupole.** p_θ is no longer conserved, and confinement over 20 µs drops
by 24–49 percentage points (paired, same initial ions).
- **Where:** almost all extra losses reach the outer wall at |z| ≈ 0.5–1.1 m,
  beyond the FRC ends.
- **Pattern:** two azimuthal sectors 180° apart in each half of the tube,
  rotated 90° between z > 0 and z < 0. This is the field-line divergence of
  a straight quadrupole superposed on B_z. See `loss_maps.png`.
  - The wall losses in each half have an m=2 amplitude |⟨e^{2iφ}⟩| of
    0.92–0.95 (noise 0.04–0.06).
  - The sector centres sit at φ ≈ 10–16° (+180°) for z > 0 and at
    −74° to −80° (+180°) for z < 0, with the bars at 45° + k·90°.
- **Size of the effect:** this is the particle-loss cost of n=2
  stabilization, and here it is large.
- **Likely an upper bound:** the model superposes the vacuum quadrupole on
  the plasma, with no screening by the FRC or the coil.

**Inner tube.** The inner separatrix sits 0.4 cm (inner tube) or 0.1 cm
(flux core) outside the 3.81 cm tube, less than one gyroradius.
- 13–23% of the ions strike the inner tube within about 1 µs of loading.
- The annular device therefore loses ions to its inner tube. The flux core
  pushes the inner separatrix closer to the tube, yet halves the prompt
  loss, presumably because its stronger field near the tube magnetizes ions
  there.

**Outputs** (outputs/particle_tracking/ions_RUN; not committed):
- `summary.json`: all numbers, including confinement curves by class, loss
  energies and times, and the m=2 loss amplitude per half;
- `confinement.png`, `loss_maps.png`, `orbit_classes.png`;
- `ions_*.npz`: per-ion initial state, loss data and class;
- `particles.json` (also written to RUN/): for the viewer.

## Multipole strength and order (`multipoles`, outputs/particle_tracking/multipoles_fine_uniform_h2)
**Question.** What does n=2 stabilization cost in ions, against multipole
strength and order? The stability threshold itself is not computed here; it
needs the 3D hybrid code (frc-ibl.11.4).

**Setup.**
- Option A, the same 4000 H⁺ (seed 1) as `ions`, 10 → 30 µs, dt = 4 ns in
  every case.
- Static multipoles (`multipole(I, order)`), bars at r = 0.20 m, 2.2 m
  long, in series through end arcs:
  - q = 4-bar quadrupole;
  - o = 8-bar octupole.
- The level is the mean |B| on the wall circle (r = 15.24 cm) at the
  midplane, as a percentage of B_e = 0.245 T.
- Changes are paired against the no-multipole case. The standard error is
  the paired one, √((b + c − (b − c)²/n)/n²).
- q7 reproduces the earlier `ions` run exactly (71.4%).

| Case | Bar current | Mean \|B\| at the separatrix / max on the wall (% B_e) | Confined at 20 / 30 µs | Change at 30 µs (points) | Lost to wall / ends | Wall-loss \|z\|, 10–50–90% | Wall-loss amplitude z>0 / z<0 (noise) | Sector centres z>0 / z<0 |
|---|---|---|---|---|---|---|---|---|
| none | – | – | 96.3 / 95.0% | – | 1.2 / 3.8% | 0.82–0.95–1.04 m | – (48 ions) | – |
| q3 | 2.34 kA | 2.2 / 4.2 | 95.3 / 87.9% | −7.0 ± 0.5 | 6.8 / 5.2% | 0.74–0.85–1.02 m | m=2 0.89 / 0.89 (0.09) | 24° / −66° |
| q5 | 3.90 kA | 3.6 / 7.0 | 93.2 / 79.5% | −15.4 ± 0.6 | 15.5 / 4.9% | 0.66–0.75–0.97 m | m=2 0.91 / 0.91 (0.06) | 21° / −71° |
| q7 | 5.46 kA | 5.1 / 9.8 | 90.9 / 71.4% | −23.6 ± 0.7 | 25.4 / 3.3% | 0.61–0.69–0.94 m | m=2 0.92 / 0.92 (0.04) | 16° / −74° |
| q9 | 7.03 kA | 6.5 / 12.6 | 88.4 / 64.0% | −30.9 ± 0.8 | 33.1 / 2.9% | 0.52–0.65–0.90 m | m=2 0.93 / 0.93 (0.04) | 14° / −76° |
| o7 | 4.83 kA | 2.9 / 7.8 | 95.9 / 90.2% | −4.8 ± 0.4 | 4.4 / 5.4% | 0.71–0.83–1.03 m | m=4 0.83 / 0.84 (0.11); m=2 0.04 / 0.03 | 16° / −31° (mod 90°) |

Sector centres are the phase of ⟨e^{imφ}⟩/m, so they repeat every 180°
(quadrupole) or 90° (octupole). Quadrupole bars sit at 45° + k·90°, octupole
bars at 22.5° + k·45°.

**Readings.**
- **Quadrupole cost grows almost linearly with strength.** It is about
  2.3 points per % of B_e at 3%, and 3.4 points per % from 5 to 9%.
  - No strength is free: 3% still loses 7 extra points by 30 µs.
  - Losses are still accelerating at 30 µs. Over 25–30 µs, q7 loses
    2.0 points/µs, q3 0.84, o7 0.55, and none 0.07. A 50–150 µs FRC would
    pay more than these 20 µs numbers.
- **The octupole is much cheaper.**
  - At the same 7% wall mean it costs 4.8 points, less than the 3%
    quadrupole.
  - Its interior field falls as r³, so at the separatrix it is 2.9% of B_e.
    Its field is also more even around the wall: max/mean is 1.11, against
    1.40 for the quadrupole.
  - At equal separatrix field it still costs less than half as much. The
    quadrupole curve gives ≈ 11 points at 2.9% at the separatrix (≈ 4% at the
    wall).
- **Loss pattern.** Both orders leave the same field-line-divergence imprint:
  2·order sectors per half, rotated 180°/order between z > 0 and z < 0
  (measured: 90° for the quadrupole, 47° for the octupole).
  - The quadrupole's m=4 amplitude (0.64–0.74) is only the harmonic of its
    narrow m=2 sectors.
  - As the quadrupole strengthens, losses move inward along the wall (median
    |z| from 0.85 to 0.65 m). Ions that would have left through the ends
    reach the wall first, so end losses fall from 5.2% to 2.9%.
- **Limitations.**
  - The multipole is a vacuum superposition: there is no screening by the
    plasma or the coil. The losses are therefore upper bounds.
  - Test particles move in fixed 2D fields: no collisions, no E_r and no
    rotation.
  - Whether any of these levels holds off n=2 is not shown. In particular,
    whether an octupole does at 2.9% at the separatrix needs frc-ibl.11.4.

**Outputs** (not committed): `summary.json` (per case: spec, confined
fraction every 5 µs, loss split, wall-loss |z| quantiles, m=2 and m=4
amplitudes, paired change), `ions_<case>.npz`, `multipoles.png`.

## RMF electrons (`rmf`, outputs/particle_tracking/rmf_even)
**Setup.**
- 300 thermal electrons loaded as above at 10 µs, in the saved fields plus
  an even-parity antenna: 1 MHz, 5 mT at the centre.
- Antenna: 120° saddles at r = 0.16 m over |z| ≤ 0.6 m, 1.12 kA per pair.
- A control run has no RMF.

**Diagnostics.**
- Each electron's θ̇ (2 ns samples, first 25% skipped) is averaged and
  compared with ω.
- Electrons are split by whether their time-averaged |B_z| is below the
  RMF amplitude.

**Approximations.** The antenna's vacuum E_z (up to ~ωB_ω r ≈ 3 kV/m) is
unscreened. In a real plasma, the electrons' co-rotation and the plasma
currents cancel it.

**Results** (2 µs = 2 RMF periods, ω dt_ce = 0.3, dt = 6.7 ps, 2.97×10⁵
steps; no electron reached a wall).

| Median θ̇/ω, by time-averaged radius | No RMF | RMF 1 MHz, 5 mT |
|---|---|---|
| 3.8–5.7 cm (22–38 electrons) | −0.06 | 0.10 |
| 5.7–7.6 cm (188–198) | 0.01 | 0.23 |
| 7.6–9.5 cm, around the null R ≈ 8.3 cm (57–84) | 0.09 | **0.93** |
| All 300 | 0.018 | 0.46 |
| Within ±30% of ω | 3% | 22% |
| Median energy at 2 µs (start: mean 129 eV) | 101 eV | 429 eV |

The paired shift of θ̇/ω (the same electrons with and without the RMF) has
median 0.31, with quartiles −0.47 and 1.11.

**Reading.**
- **Entrainment at the null.** Electrons near the null co-rotate with the RMF
  at essentially ω, in the +θ direction that carries the FRC's diamagnetic
  current. This is the single-particle mechanism of RMF current drive:
  where the RMF dominates B_z, the E×B drift in the rotating field is the
  perpendicular part of rigid rotation at ω.
- **Strong B_z elsewhere.** Inside, where B_z ≫ B_ω, test electrons cannot
  co-rotate. That would need the radial E_r of a co-rotating plasma, which
  test particles do not have. Instead they are driven along B by the
  unscreened E_z and heated about fourfold in 2 µs.
- **What a test particle cannot show.** Full RMF penetration and sustained
  current are self-consistent plasma effects, left to the hybrid model
  (frc-ibl.11.5).

**Frequency ordering.** At 1 MHz the ordering ω_ci < ω < ω_ce holds only
where |B| < 0.066 T:
- f_ci at B_e = 3.7 MHz and f_ce = 6.9 GHz;
- ω > ω_ci at 47% of the loaded electron positions, around the null.

At 0.2 T an H⁺ RMF must run above ~4 MHz (PFRC-like, odd parity possible)
to keep the ions from following it everywhere. The TCS range of 85–240 kHz
belongs to ~0.01–0.02 T devices.

**Control-run artefact.** A few unmagnetized electrons at the null run away
in the resistive E_θ, which is unopposed because there are no collisions.
They lift the control's mean energy to 228 eV, while its median stays at
101 eV.

## Viewer trails
`ions` writes RUN/particles.json, which `formation_display.py` copies into
its output. It holds actual integrated orbits, positions every 20 ns:
- up to five confined ions of each class;
- five lost ions.

It has one set per field case (no quadrupole, quadrupole at the wall, at
the separatrix).
