# Virtual plasma bench: project understanding and present status

The goal is a scriptable, eventually interactive plasma laboratory: define a
chamber, gas/plasma inventory, coils and prescribed waveforms; evolve the plasma;
inspect fields, currents, temperatures, motion, energy and magnetic topology.
FRC formation is the first major application. Field reversal, closed flux,
force balance, lifetime, three-dimensional stability and particle confinement
are separate claims with separate evidence requirements.

The user's provisional envelope is about 1 m long, 15 cm diameter and
0.1–100 mTorr. Available gases include argon, H2, D2, N2 and O2. These are inputs
to investigate, not established operating settings. Coil hardware remains an
abstract prescribed waveform; this repository does not design pulsed power.

## What exists

| Area | Implemented work | Present boundary |
|---|---|---|
| Vacuum fields and drive | Axisymmetric circular coils, prescribed waveforms, flux and field calculations, magnetic-only exterior coupling | Actual exterior vacuum is retained; chamber/end/contact assumptions vary by model and must be read with the chosen experiment. |
| Moving plasma | Conservative material motion, separate electron/ion temperatures, finite carrier masses, Hall/current inertia, full axisymmetric field/current components and moving magnetic vacuum | Several research paths exist. They are not one finished, validated formation solver. Current, thermal and edge accuracy remain open. |
| Configured coil experiments | Metre-scale JSON experiments and complete bias/pulse trajectories, including paired finite-drift Coulomb force/heat | Preparation is assumed. Strong-pulse current is not time/spatially converged; a refined trajectory stopped at the unimplemented wall/sheath contact. |
| Atomic physics | Attributed argon ionization/recombination, resolved-level source studies, radiation components, conservative momentum/energy sources and sourced electron/ground-neutral elastic moments | Product branching, other neutral channels, photon transport and complete driven coupling remain absent. |
| Changing species in motion | Separate-support planar Euler–Poisson transport plus actual argon reactions, exact empty exterior and explicit energy accounting | Newly corrected runs finish, but electron-current refinement and local thermal transfer remain inaccurate. This planar component is separate from the coil experiment. |
| Molecular gases | New sourced H2/N2/O2 total-ionization cross sections and bounded Maxwellian integrals | Total ionization is only one component. D2 data, products, internal energy, dissociation, attachment and complete transport/coupling remain unfinished. |
| Visualization | Scientific plots, topology diagnostics and interactive 3D playback of saved axisymmetric runs | Playback shows its archived run. It is not evidence of a three-dimensional solver or an interactive, validated lab-design tool. |

See [the configured experiment](MATERIAL_EXPERIMENT.md),
[finite-drift collision results](FINITE_DRIFT_COLLISIONS.md),
[atomic source](ATOMIC_NETWORK.md) and [gas scope](MOLECULAR_GASES.md).
The earlier H2 inventory option assumes complete dissociation into fully ionized
protons; it is not a molecular hydrogen discharge model.

## What the latest continuation established

The active goal continuation adds conservative variance transfer and explicit
mechanical refinement. Refining the actual trace-ion transport stage reduces
its peak temperature from7.04 to0.35eV while retaining every population and
material cut. Complete20ns enriched10ns-source runs finish, but the shorter
source-interval case exposed arithmetic underflow in a cold trace-cell sound
calculation and loss of a small local mixing-heat contribution. Equivalent
ratios and explicit binary number/energy scaling address these arithmetic
errors; full shorter-interval evaluation remains active. Current
still differs substantially from an existing fine numerical trajectory, so a
completed128-electron/32-heavy run reproduces the earlier128-cells-per-species current
within0.056% at lower observed cost. It remains34.2% away from the existing
1024-cell numerical reference. The completed256-electron/32-heavy run reduces
that current difference to19.1% and electron-front difference3.39→1.86mm.
The1024reference is itself unconverged; current accuracy is still open. Independent analytic
expansion families now complete their original full durations. [VARIANCE_TRANSFER.md](VARIANCE_TRANSFER.md) records all
current evidence and limitations. Defaults remain unchanged.

The original variance5ns trajectory now also completes20ns. Its10/5ns current
difference is0.0434%, but Ar8 T differs0.1612eV RMS, keeping local accuracy open.
Charged collisions now feed back into the moving reaction slab in a complete
20ns run. Sourced electron/ground-neutral elastic moments complete independent
full relaxation comparisons, including recoil and heat, and their moving source
completes20ns with total-energy drift2.842e-11. See [neutral-collision evidence](ELASTIC_NEUTRALS.md)
and [moving charged-collision evidence](ATOMIC_COLLISIONS.md). These are distinct
component advances, not a completed laboratory simulator.

The full3.5us heated source now includes electron/ground-neutral collisions at
two time resolutions. It shows why vector current and neutral motion matter:
azimuthal current damps faster while different neutral/ion flows sustain axial
current. Ground-ion/neutral scattering is the next missing physical channel;
the published table's low-energy gap covers up to11.6% of the modeled pair
distribution. Excited-neutral channels are also incomplete. The coupled
source's small time-halving error does not remove these physical limitations.
The first uniform-density electron experiment completes but accidentally selects
only102velocity cells through the old cosine targets. Its worsened current/front
comparison rejects that mixed mesh. The constructor call is corrected, and the
new full-duration run starts with128 actual uniform velocity cells; its accuracy
comparison remains pending. Defaults retain the cosine mesh.

The extreme trace-ion velocity was an approximation error across density
contrast, reproduced independently at 100-digit precision. A constrained
velocity transfer fixes that specific failure while preserving block momentum
and each fine cell's total energy. It does not establish local thermal accuracy.
Affine edge extension also removes artificial tiny velocity elements; independent
full neutral/charged expansion families complete. Exact force integration and
coordinate-search improvements reduce cost without changing the physical model.

Complete corrected 32/64/128-cell argon trajectories reach the same 20 ns target.
Successive electric-field differences decrease from 6.65% to 2.70%, while current
differences increase from 32.21% to 38.94%. Full controls with reactions and
remaps disabled reproduce these differences. At 64/128, about 87.7% of the
electron-current squared difference lies inside |x|<9 mm and the velocity term
dominates its exact density/velocity decomposition. This points to the existing
transport approximation, rather than reaction return, as the dominant spatial
current limitation in this case. Conservation alone does not establish accuracy.

The corrected 5 ns-reaction-interval trajectory also completes the full 20 ns,
retaining 2.42 million density cells. Relative to 10 ns reaction intervals,
E/J differences are 0.000921% / 0.0379%. The newborn velocity-mesh defect is
resolved, but Ar8+ temperature differs by 0.22056 eV RMS and 3.34554 eV locally.
The complete result therefore keeps thermal accuracy and the broader goal open.

The H2/N2/O2 component reproduces all 293 NIST theoretical table points within
printed rounding and agrees with independent rate quadrature. That verifies the
implementation of this component; it is not experimental discharge validation.
Detailed measurements, failures and source archives are linked from
[reaction-transfer results](REACTION_TRANSFER.md) and [gas scope](MOLECULAR_GASES.md).

The 20 ns planar case uses 1e14/m³ initial component densities and a 20 mm slab.
It is a numerical component problem, not the user's neutral-fill pressure,
full chamber geometry or demonstrated preparation procedure.

## What remains before lab-facing predictions

1. Resolve local current/front and thermal accuracy in the moving species solver,
   including complete source-time sensitivity with every species retained.
2. Couple the physical reaction, collision and heat models into the driven
   moving-plasma path, with gas-specific chemistry and distribution-validity limits.
3. Implement the required chamber/wall/sheath/end interactions and preparation
   physics, rather than continuing through an unsupported contact or assuming it.
4. Complete formation trajectories with converged fields/topology and compare
   matching physical inputs against independent experimental measurements.
5. Expose the supported model and diagnostics through the interactive interface
   without presenting axisymmetric results as proof of 3D stability/confinement.

The project is substantial working research software, but is not finished or
validated for choosing laboratory construction parameters. The full goal stays
open. Run status, exact next action and local checkpoint state are in
[THREAD_HANDOFF.md](THREAD_HANDOFF.md); older dated entries are history. Work is
saved locally, without a configured remote or publication. Verification follows
complete simulations and independent comparisons, not the historical unit suite.
