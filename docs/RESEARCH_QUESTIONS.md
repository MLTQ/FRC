# Research queue

- [ ] Choose a formation-time conductivity closure; test anisotropy and Hall terms.
- [ ] Connect gas pressure, ionization fraction, density and pre-ionization histories.
- [ ] Validate electron-neutral and ion-neutral cross sections for argon and other gases.
- [ ] Establish collisionality / kinetic limits with defensible mean free paths.
- [ ] Identify where di/L, gyro scales and pulse times require Hall/two-fluid physics.
- [ ] Replace distant fixed induced-flux boundary with a verified open-boundary map.
- [ ] Quantify end effects, actual quartz geometry and any conductive structures.
- [ ] Refine interface placement, temporal error and spatial convergence near skin layers.
- [ ] Verify topology with manufactured islands, axis nulls, separatrices and noise.
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
- [ ] Replace the fixed chamber interface with a defensible evolving free-plasma
  treatment without exterior fluid floors; validate gas loss and actual end conditions.
- [ ] Add thermodynamic feedback before interpreting large Ohmic heating as plausible.
- [ ] Test nonaxisymmetric stability only after a 3-D solver exists.

First comparisons use equal peak absolute vacuum perturbation at the chamber
center and identical geometry, bias, conductivity and output window. They are
not equal source energy, equal Ohmic heating, or equal positive field excursion.
All example dimensions, densities, currents and conductivities are exploratory
inputs, not inferred laboratory requirements or hardware recommendations.
