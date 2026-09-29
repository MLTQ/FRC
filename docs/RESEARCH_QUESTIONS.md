# Research queue

- [ ] Choose a formation-time conductivity closure; test anisotropy and Hall terms.
- [ ] Connect gas pressure, ionization fraction, density and pre-ionization histories.
- [ ] Validate electron-neutral and ion-neutral cross sections for argon and other gases.
- [ ] Establish collisionality / kinetic limits with defensible mean free paths.
- [ ] Identify where di/L, gyro scales and pulse times require Hall/two-fluid physics.
- [x] Add a verified open magnetic boundary using decaying multipoles, including
  full-space magnetic energy, asymmetric source checks and matching-surface studies.
- [ ] Verify and, if necessary, improve open-boundary efficiency for long chambers;
  the current single-center expansion requires a source-enclosing sphere.
- [ ] Quantify end effects, actual quartz geometry and any conductive structures.
- [ ] Refine interface placement, temporal error and spatial convergence near skin layers.
- [ ] Verify topology with manufactured islands, axis nulls, separatrices and noise.
  A clean analytic spherical FRC now verifies the existing O-ring/axis-null
  detector and approximately second-order separatrix-volume convergence.
  Noisy, degenerate and subgrid cases remain open.
- [ ] Define an FRC-like transient using independent reversal, topology and force-balance criteria.
- [ ] Define excluded-flux radius against a specified probe geometry/reference flux.
- [ ] Define equal-energy waveform comparisons including source/plasma mutual energy.
- [ ] Add physical B-dot, flux-loop, Rogowski and optical diagnostic response models.
- [x] Implement a conservative bounded-domain resistive-MHD reference kernel;
  verify cylindrical fast waves, 2-D advection, shock propagation and heating.
- [x] Reduce smooth-field numerical diffusion with limited second-order reconstruction;
  verify annular moments, axis forces and the full one-dimensional UCT limit.
- [x] Add an independent curved closed-flux Grad–Shafranov equilibrium and quantify drift.
- [x] Add consistent positive face/corner reconstruction for low-pressure states;
  retain disabled-reference failures, unchanged cell budgets and limiter counts.
- [ ] Improve cold rarefaction pressure accuracy and test stronger full-update
  robustness; positive traces alone are not a proof of accurate evolution.
- [ ] Improve curved-equilibrium force balance and quantify residual artificial flow.
- [ ] Couple dynamic MHD to prescribed coils, quartz boundaries and genuine vacuum.
- [x] Verify an intermediate fixed, gas-filled chamber with moving fluid,
  prescribed coils, magnetic-only exterior and a discrete source-work budget.
- [x] Separate elliptic solver residual from finite magnetic-domain error using
  an independent infinite-space circular-loop integral and exact-boundary checks.
- [ ] Resolve pulse-driven skin layers and compression in the coupled reference;
  source work is nearly converged before peak flow is converged on current meshes.
- [ ] Establish magnetic-shock/strong-compression accuracy for the split coupling;
  energy pairing and positive internal energy are insufficient evidence.
- [x] Expose centered magnetic-shock oscillations with an independent cylindrical
  Rankine–Hugoniot/steady-flow reference in the actual two-dimensional kernels.
- [x] Add selectable vanishing compressive numerical magnetic diffusion with
  separately measured heat; verify the faster cached-vacuum implicit solve
  against the retained direct method under both magnetic boundary choices.
- [ ] Resolve the remaining irregular B/rho error and overshoot, establish local
  conservative shock limits, and extend shock tests to other strengths and
  orientations. Numerical diffusion is not physical transport or reconnection.
- [x] Add unsplit conservative fluid/UCT coupling to the existing magnetic vacuum;
  verify original shocks, independent smooth solutions and driven source-work
  refinement, while distinguishing fluid conservation from full-space accuracy.
- [x] Add a declared whole-update conservative flux blend when positive traces
  alone fail, preserving budgets and CT with no cell floor or energy reset.
- [x] Localize complete-update fallback by shared face/corner flux replacement,
  rechecking newly invalid neighbors. Original nr=192 shock B error improves
  about 20%, with <0.86% cells touched in any stage and none in the central
  reference window; no cell floors or energy repairs.
- [ ] Resolve the independent whole-space shock energy error exposed by local
  limiting. Central profile accuracy improves but nr=192 response error grows
  to 1.40e-4 of initial energy, then 5.72e-4 at nr=384 (30.20% of net outward
  fluid energy); exact fluid conservation is insufficient.
  Check finer meshes, time step, interface flux consistency and boundary effects.
- [x] Derive and independently perturbation-check the instantaneous response
  energy rate at common states, separating boundary fluid flux, cell magnetic
  energy and full-space nodal work. The audit does not alter evolved diagnostics.
- [x] Add compatible staggered fluid magnetic energy and natural boundary
  Poynting work; independently verify the constant-bias semidiscrete identity.
- [x] Separate remaining finite-step error with a quadratic RK stage audit and
  verify optional SSP RK3 using independent nonautonomous and actual-MHD controls.
- [ ] Validate the experimental fixed first-order boundary strip under source
  pulses, skin-layer refinement, matching/end changes and finer shock meshes.
  It preserves the sharper second-order interior in the initial controlled probes;
  its explicit order reduction is not yet a production default or physical closure.
- [ ] Refine and validate magnetic interface traction and numerical heating;
  conservative pulse response error is 0.1021% of source work at nr=64, but
  numerical interface heat and driven peak flow are still resolution-sensitive.
- [ ] Replace the fixed chamber interface with a defensible evolving free-plasma
  treatment without exterior fluid floors; validate gas loss and actual end conditions.
  The independently verified spherical Hill-vortex reference supplies a true
  vacuum, continuous boundary traction and exact full-space energy target;
  its numerical evolution and a thermodynamic closure are not yet implemented.
- [ ] Add thermodynamic feedback before interpreting large Ohmic heating as plausible.
- [ ] Test nonaxisymmetric stability only after a 3-D solver exists.

First comparisons use equal peak absolute vacuum perturbation at the chamber
center and identical geometry, bias, conductivity and output window. They are
not equal source energy, equal Ohmic heating, or equal positive field excursion.
All example dimensions, densities, currents and conductivities are exploratory
inputs, not inferred laboratory requirements or hardware recommendations.
