# Research queue

## User experiment envelope (2026-10-04)

The user proposes a long quartz tube, approximately **1 m long and 10–20 cm
diameter**, with dimensions open to revision. Gas species is selectable; the
user has several gases and pumping options. The proposed **neutral fill
pressure is 0.1–100 mTorr** (about 0.0133–13.3 Pa), not an electron or ion
number density. Gas temperature and ionization fraction remain independent
inputs; for an ideal gas n_neutral=p/(k_B*T_gas), while carrier density depends
on ionization and charge state. Existing smaller argon verification fixtures
are not recommendations for this experiment.

The bias field will be built for the experiment; a Helmholtz-like coil
arrangement is preferred if suitable. No field magnitude is fixed. Compare
prescribed coil fields/waveforms and their spatial uniformity across the
selected plasma region before selecting a field range. The 1 m tube length
does not imply a 1 m uniformly magnetized or fully ionized plasma column.
Do not infer a validated gas, pressure, pulse or hardware design from this
envelope: model validity, attainable reversal/topology and numerical convergence
must determine the useful region of parameter space. Hardware stays abstract.

## Open and completed research

- [ ] Choose a formation-time conductivity closure; test anisotropy and Hall terms.
- [x] Add independently checked conditional Coulomb coefficients and local
  Hall/electron-inertia/Debye, drift and collisionality diagnostics; retain
  explicit assumptions and never turn an absent warning into a validity claim.
- [x] Implement a fixed-material current-inertia component with physical electron
  mass, variable carrier/resistivity profiles and true magnetic vacuum. Verify
  continuum decay/screening, variable coefficients, current kinetic energy,
  positive heat and independent time order; retain coarse false-reversal evidence.
- [x] Implement and verify a controlled axisymmetric linear Hall/toroidal
  reference, both helicities, physical electron inertia, damping and independent
  cylindrical/parallel dispersion (49 cylinder cases, 16 parallel modes).
- [x] Integrate conservative toroidal field/angular momentum into a nonlinear
  bounded compressible core, with independent stresses, finite-amplitude
  transport, waves and compatible scalar resistive work (33 complete cases).
- [x] Derive and independently verify finite-sound-speed Hall modes and a
  finite-amplitude variable-density Hall/Poynting snapshot; record the complete
  42-case experimental spatial study with global work and local errors.
- [x] Correct the experimental Hall operator's nonvanishing first-annulus gdot
  error (21.46% through nr128) with consistent regular-axis energy/current/
  electric moments, then verify local transport and stable actual evolution.
- [x] Couple nonlinear Hall to the compressible dynamic state with explicit
  compatible magnetic energy, nonlinear midpoint and physical resistive heat;
  verify three actual wave branches, time/mesh controls, nonlinear flow and
  independent conducting Bessel decay (63 complete cases, 257 passing tests).
- [ ] Improve physical magnetic representation accuracy: nr32 wave perturbation
  errors remain about 4%, dominated by the initial lowest-order field space.
- [x] Derive and verify conservative dual charge flux, electron pressure work
  and entropy/enthalpy drift as an isolated component; retain full wall trace
  and discrete physical entropy alongside energy. All 46 records complete,
  including actual heat translation, nonlinear battery and a strict failure.
- [ ] Couple this dual electron thermal state to moving primal ion mass and
  total energy with a consistent boundary model; then warm waves, collision
  exchange and full finite-mass current terms. TWO_TEMPERATURE_PLAN retains the
  conservative mass-map and local energy-flux requirements. Isolated electron
  pressure dynamics do not supply a physical heated edge or formation model.
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
- [x] Measure the experimental fixed first-order boundary strip under source
  pulses, analytic skin-layer refinement, matching/end changes and nr=384 shocks.
  It preserves interior shock accuracy and reduces fine shock response error
  to 1.01130e-7 J, but changes nr=128 driven velocity by up to 9.767% versus
  the standard method. It remains experimental and is not a physical closure.
- [ ] Improve boundary/flow accuracy using the new full-field evidence;
  64→128 velocity changes are 7.059% standard and 10.423% with the strip.
  Small source-budget error and converged matching/time controls do not remove
  this spatial limitation. See BOUNDARY_VERIFICATION.md.
- [ ] Refine and validate magnetic interface traction and numerical heating;
  compatible nr=128 pulse response error is 0.000838% of source work, but
  flow, density perturbations and current remain resolution-sensitive.
  Local numerical eta and integrated numerical/physical heat are separately
  reported; the latter excludes implicit Rusanov/UCT diffusion.
- [ ] Replace the fixed chamber interface with a defensible evolving free-plasma
  treatment without exterior fluid floors; validate gas loss and actual end conditions.
  The independently verified spherical Hill-vortex reference supplies a true
  vacuum, continuous boundary traction and exact full-space energy target;
  its ideal material-mesh evolution now has measured force and drift refinement.
  Scalar resistive coupling is now verified separately; physical surface
  thermodynamics and resolved nonideal moving flow remain open.
- [x] Implement a two-dimensional ideal material mesh with gas only inside an
  evolving surface, regular a=psi/r² magnetic energy, analytic geometry forces,
  an exactly curved matching sphere and full-space decaying vacuum modes.
- [x] Verify actual free gas expansion, initialized FRC force/drift, independent
  sphere energy/virtual force, prescribed uniform-bias work, and material flux
  interpolation error. See MOVING_PLASMA_VERIFICATION.md for recorded studies.
- [x] Extend the moving reference to resistive induction and positive physical
  Ohmic heating with compatible full-space energy and entropy accounting;
  independently verify diffusion/current/heat, source work and time splitting.
- [x] Expose raw midpoint false reversal in an independently diagonalized decay
  mode; add adaptive field/current/local-heat controls with positive physical
  heating and atomic failure. Local tolerance is not a global topology guarantee.
- [ ] Resolve moving resistive flow and surface thermal behavior: 16→32 material
  grids differ by 63.821% in driven velocity at D=0.1 m²/s, despite much smaller
  time/boundary controls and about second-order fixed-sphere current/heat errors.
  D=1 changes the physics and still differs by 53.915%; it is not a numerical
  remedy. Inspect spatial modes, axis/interface forces and geometry before
  interpreting initialized reversal or improved global energy accounting.
- [x] Add explicit radial grading on both sides of the material surface and
  independent angular/time/tolerance/matching controls; retain the unchanged
  uniform default and measure interior-resolution tradeoffs with exact gas
  expansion and sphere-decay references.
- [x] Compare differently graded states on common initial physical material
  points, integrating exact cell intersections rather than missing discontinuous
  cell edges. Check independent clipped moments/volume and quadrature order.
- [ ] Resolve the conditional thermal compatibility defect in the heated sharp
  interface. For continuous B and a smooth material p=0 trace with bounded
  div(v), the model requires eta J²=0 at the surface; the initialized ideal FRC
  has nonzero J there and violates this at the first time derivative. This is
  not a proof about all weak solutions. Retain the original stress case, derive
  physically justified edge/current/thermal assumptions, and verify the required
  non-smooth dynamics or extended transport model. See SURFACE_COMPATIBILITY.md.
- [ ] Couple spatial circular-coil sources to moving material coordinates,
  including geometric source derivatives and source-work controls; the new
  material solver currently supports only prescribed uniform axial bias.
- [ ] Address material-mesh distortion/low-energy modes, shock entropy, remapping
  and continuum frozen flux between vertices. Preserved nodal psi alone is not
  a topology guarantee; nonuniform-deformation flux interpolation error is
  measured explicitly and must remain separate from physical resistivity.
- [ ] Reduce and quantify the finite magnetic mesh's spurious axial force on a
  translated plasma and its matching-surface dependence under stronger motion.
- [ ] Add thermodynamic feedback before interpreting large Ohmic heating as plausible.
- [ ] Test nonaxisymmetric stability only after a 3-D solver exists.

First comparisons use equal peak absolute vacuum perturbation at the chamber
center and identical geometry, bias, conductivity and output window. They are
not equal source energy, equal Ohmic heating, or equal positive field excursion.
All example dimensions, densities, currents and conductivities are exploratory
inputs, not inferred laboratory requirements or hardware recommendations.
