# 3D hybrid model of the formed FRC (frc-ibl.11)

## Why
The 2D formation solver cannot show the effects that set FRC lifetime or the
user's experiments around it:
- n=2 rotation and tilt;
- quadrupole stabilization;
- rotating-magnetic-field (RMF) current drive;
- particle confinement.

FRC ions in the 12" device are kinetic: s̄ ≈ 2–3, only a few gyroradii across
the separatrix. The standard self-consistent model for this regime is
**hybrid**: kinetic ions with fluid electrons. Formation stays in the
benchmarked 2D solver; the 3D model starts from a formed FRC.

## Model
**Ions.** Full-orbit macro-particles (H⁺ or D⁺), advanced by the Boris push
with trilinear (CIC) weighting. Density n and flux n u_i are deposited on the
grid nodes.

**Electrons.** A massless, quasi-neutral fluid (n_e = n_i) obeying the
generalized Ohm's law, with the Hall term and electron pressure:

    E = −u_i×B + J×B/(n e) − ∇p_e/(n e) + η J − η_H ∇²J,   J = ∇×B/μ0

- Electron pressure: isothermal (p_e = n T_e) or adiabatic.
- η is physical resistivity (uniform, or Spitzer later).
- η_H is a **declared** hyper-resistivity. It damps grid-scale whistlers and
  grid modes, and its Joule heat is reported.

**Fields.** Darwin, with no displacement current.
- Faraday's law, ∂B/∂t = −∇×E, on a Yee grid: B on faces, E and J on edges.
- ∇·B is preserved to round-off.

**Low density.** Ohm's law uses max(n, n_floor). The floor is declared and its
use is counted. A vacuum region is not a plasma; see the geometry milestone.

## Time integration
Per ion step dt:
1. E^n is computed from Ohm's law with B^n, n^n and u^n. Here
   u^n = 1.5 u^{n+½−1} − 0.5 u^{n−½−1} is extrapolated from the two most
   recent half-step fluxes.
2. Boris push: v^{n−½} → v^{n+½} with E^n and B^n at x^n, then
   x^{n+1} = x^n + dt v^{n+½}.
3. Deposit n^{n+1} at x^{n+1}, and the flux at x^{n+½} with v^{n+½}.
4. Advance B from t^n to t^{n+1} in subcycles of SSP-RK3. Ohm's law is
   re-evaluated every stage with n interpolated in time and u = u^{n+½}.
   RK3's stability on the imaginary axis (|λ dt| ≤ √3) covers the dispersive
   whistler.

**Time-step limits:**
- ions: dt ≤ 0.1/Ω_i, and a particle travels < 0.5 cell per step;
- fields: dt_B ≤ √3 / ω_max, with ω_max = v_A k_max (1 + k_max d_i) at the
  grid scale.

## Geometry (milestones)
- **M1:** periodic box, for verification.
- **M2:** the 12" device, on a Cartesian grid, which avoids an axis
  singularity.
  - The glass tube (r_w = 15.24 cm) and its ends absorb ions.
  - The single-turn coil is a perfect conductor at r_c = 16.5 cm, with
    conducting planes beyond the tube ends, as in the 2D model.
  - Exact vacuum fills the gap between tube and coil (no ions; J = 0). How it
    is enforced is decided in M2 and declared. Any approximation is measured
    against the 2D solver's exact vacuum.
- **M3:** start from the formed FRC of the 2D solver (flux, density, pressure,
  velocity) and from an analytic rigid-rotor equilibrium. Check that the n=0
  equilibrium holds, and compare n=0 evolution with the 2D run.
- **M4:** n=2 rotation and tilt; a quadrupole at 7% of B_e.
- **M5:** RMF antennas, even and odd parity.
- **M6:** hybrid runs and ion trails in the 3D viewer.

## Verification (independent, no unit tests)
- **Particle orbits:** the Boris gyration phase error against the exact orbit,
  and exact energy in static fields.
- **Parallel waves** in a uniform cold plasma, against the exact Hall-MHD
  (cold-ion hybrid) dispersion:

      ω² = k² v_A² (1 + ω/Ω_i)

  This gives the whistler branch (ω > 0) and the ion-cyclotron branch
  (ω → −Ω_i). Measure the frequency for k d_i = 0.25–2 and its convergence
  with the grid.
- **Conservation:**
  - ∇·B to round-off;
  - total energy (ion kinetic + magnetic, cold electrons, η = 0) against step
    size.
- **Later:**
  - the rigid-rotor equilibrium (M3);
  - published n=2 and tilt growth rates (M4);
  - RMF penetration and current against theory (M5).

## Performance
Particle push and deposit are JIT-compiled (numba, parallel over 10 cores);
the field operators use numpy. Target scale: about 64×64×256 cells with 10–50
particles per cell for the 12" device.

## Code
- `plasma_bench/hybrid3d.py`: grid, Yee operators, Ohm's law, Boris push,
  deposit, stepper.
- `examples/hybrid3d_verification.py`: the checks above.
