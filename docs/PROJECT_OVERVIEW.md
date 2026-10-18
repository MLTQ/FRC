# Virtual plasma bench: project understanding and present status

The goal is a scriptable, eventually interactive plasma laboratory: define a
chamber, gas/plasma inventory, coils and prescribed waveforms; evolve the plasma;
inspect fields, currents, temperatures, motion, energy and magnetic topology.
FRC formation is the first major application. Field reversal, closed flux,
force balance, lifetime, three-dimensional stability and particle confinement
are separate claims with separate evidence requirements.

The user's lab envelope (2026-10-18) is a Pyrex tube about 15 cm in diameter
and 1 m long, with H2/D2 as working gases and N2 for commissioning. The goal is
now a specific buildable benchtop FRC experiment. The repository specifies the
required field and current waveforms and the equivalent coil/bank numbers; the
user builds the hardware. The current status and plan are in
[GOAL_STATUS.md](GOAL_STATUS.md); the first-cut regime is in
[FRC_DESIGN_WINDOW.md](FRC_DESIGN_WINDOW.md).

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


## Paused component work (argon)

The 1D argon reaction/transport slab and the Vlasov slab are complete,
energy-conserving component studies. Their current and local thermal accuracy
are unresolved, and they are paused under the benchtop refocus. Evidence and
limitations: [REACTION_TRANSFER.md](REACTION_TRANSFER.md),
[DENSITY_INTERFACE_TRANSPORT.md](DENSITY_INTERFACE_TRANSPORT.md),
[VARIANCE_TRANSFER.md](VARIANCE_TRANSFER.md), [VLASOV_SLAB.md](VLASOV_SLAB.md),
[ELASTIC_NEUTRALS.md](ELASTIC_NEUTRALS.md) and
[ATOMIC_COLLISIONS.md](ATOMIC_COLLISIONS.md). The H2/N2/O2 total-ionization
component reproduces the 293 NIST theoretical table points; that verifies the
implementation, not discharge behaviour.

## What remains before lab-facing predictions

1. Reproduce a documented small field-reversed theta-pinch experiment with the
   most complete coil-driven solver path (beads frc-ibl.2).
2. Add only what that benchmark requires: hydrogenic plasma, pre-ionization and
   bias-flux trapping, wall/end contact, and a coil circuit with crowbar.
3. Scan the design for the user's geometry, with converged fields and topology
   and an independent cross-check.
4. Do not present axisymmetric results as proof of 3D stability or confinement.

Verification follows complete simulations and independent comparisons, not the
historical unit suite. Work is local, with no remote.
