# particle_tracking.py

## Intent
Full-orbit test particles (H⁺, D⁺, e⁻) in the fields of a saved 2D formation
run, with analytic 3D coil fields superposed. This is step 1 toward the 3D
model (frc-ibl.12). It shows orbit classes, confinement, and wall and end
losses, including the non-axisymmetric quadrupole and RMF fields that the 2D
solver cannot represent.

Test particles do not act back on the fields. There is no self-consistent
E_r, no plasma screening of coil fields, and no collisions. The
self-consistent kinetic-ion model is the hybrid code (frc-ibl.11).

## Fields
**`FluxField(r, z, times, psi)`.** Axisymmetric fields from the poloidal flux
per radian ψ(r, z) saved at `times` on a uniform grid with r₀ = 0.
- **Representation:** ψ = r² g(r, z, t). g is a cubic B-spline interpolant,
  with mirror-even ends: even in r at the axis, which keeps g smooth there.
  The axis row uses g(0) = (16ψ₁ − ψ₂)/(12 dr²), exact for ψ = a r² + b r⁴.
- **Time:** linear in t between saved states. dg/dt is piecewise constant,
  and the spline coefficients of g + i dg/dt are blended once per call.
- **Fields:**
  - B = (−x g_z, −y g_z, 2g + r g_r): exactly ∇ × (A_θ θ̂) with A_θ = r g.
    It is divergence-free and regular on the axis.
  - E = ∂_t g (y, −x, 0), i.e. E_θ = −r ∂_t g = −(∂ψ/∂t)/r. This is exactly
    Faraday's law for the interpolated ψ.
- **Exact invariant:** p_θ = γ m r v_θ + q ψ is an exact invariant of these
  axisymmetric fields, for any time dependence of ψ. The only p_θ error is
  the integrator's.
- **Interpolation:** reproduces the saved node values to ~1e-16.

**`Filaments(paths, current)`.** Closed polylines of straight filaments
carrying I(t), where `current(t)` returns (I, dI/dt).
- **B:** finite-segment Biot–Savart,
  (μ0 I/4π) (R₁×R₂)(|R₁|+|R₂|)/(|R₁||R₂|(|R₁||R₂| + R₁·R₂)).
- **E = −(dI/dt) A:** the segment vector potential is
  (μ0 I/4π) û ln((|R₁|+|R₂|+L)/(|R₁|+|R₂|−L)). Because the paths are
  closed, A is divergence-free, giving a charge-free inductive E. There is no
  retardation and no capacitive antenna field.

**Coil builders.**
- `quadrupole(I)`: the straight quadrupole of BUILD_DOSSIER.md §6.
  - Four axial bars at r = 0.20 m and azimuths 45° + k·90°, running
    z = ±1.1 m, with alternating current directions.
  - The bars are in series, joined by 90° arcs (six segments each) at
    r = 0.20 m on alternate ends. The current path is therefore one closed
    loop that never enters the tube.
- `rmf_antenna(f, B_ω, parity)`: two orthogonal saddle pairs driven with
  I₀ cos ωt (field along x) and I₀ sin ωt (along y), giving a transverse
  field that rotates in +θ.
  - Geometry: 120° saddles at `radius` (default 0.16 m), over
    |z| ≤ `half_length` (0.6 m).
  - Odd parity splits each pair at z = 0 with opposite currents, so B_⊥
    reverses across the midplane.
  - I₀ gives |B_⊥| = B_ω at the centre (even) or at z = half_length/2 (odd).

**`Superposition(*parts)`** sums the parts. Its `psi()` comes from the
axisymmetric part.

**Verification fields.** `UniformField` and `GradientField`
(B = B0(1 + x/L) ẑ).

## Integrator (`track`, `boris_kick`)
**Momentum variable.** u = γv per unit mass.

**Boris drift–kick–drift.** Each step:
1. x ← x + (dt/2) u/γ;
2. fields at (x, t + dt/2);
3. half electric kick;
4. magnetic rotation with γ at the half-kicked momentum;
5. half electric kick;
6. x ← x + (dt/2) u/γ.

**Properties.**
- Positions and momenta are synchronous, so energy and p_θ are evaluated
  without interpolation.
- A pure magnetic kick rotates u exactly, so the energy is conserved to
  round-off in any static B.
- The gyro-phase error is the Boris angle 2 atan(ω dt/2) per step. The DKD
  orbit radius in uniform B is exactly u/(γω).
- `energy_ev` uses m u²/(γ + 1). The equivalent (γ − 1) m c² loses about 9
  digits for 100 eV ions.

**Time step.** Callers choose dt with an integer number of steps per saved
interval (examples/particle_tracking.py `time_step`). No kick then straddles
a jump in ∂ψ/∂t, and second order is kept in time-dependent fields.

## Walls and losses (`Walls`, `track`)
**Absorbing surfaces.** r = `r_outer`, r = `r_inner` (an inner tube, when
> 0) and |z| = `z_end`.

**Removal.** A particle whose step ends outside is removed from the active
arrays. Its id records:
- the crossing point and time, interpolated within the step (a quadratic for
  the cylinders, linear for the end planes);
- the surface code (`SURFACES`: wall, inner tube, end +z, end −z);
- its kinetic energy.

`observe(t, population, B_half)` is called after every step with the active
particles.

## Verification (examples/particle_tracking.py `verify`; numbers in examples/particle_tracking.md)
The pusher and fields are checked against independent analytic results:
- gyration frequency and radius;
- E×B and ∇B drift speeds;
- the mirror loss cone, from the exact flux of two current loops passed
  through `FluxField`;
- the quadrupole against four infinite line currents;
- RMF Faraday's law by central differences;
- energy conservation (frozen FRC plus the 3D quadrupole);
- p_θ conservation, with convergence order, in frozen and time-dependent
  saved FRC fields.

## Contracts
- **Units:** SI, with energies in eV. Positions and momenta are (N, 3)
  Cartesian arrays.
- `Population.ids` index the per-particle loss arrays, which keep their
  original size.
- `FluxField` raises when it is asked for a time outside its saved states.
