# Moving adiabatic electron/ion reference

`AdiabaticSpeciesMHD` extends `ToroidalMHD`'s ideal SSPRK2 with conservative
electron U advection using each stage's actual mass flux through `DualMass`.
Electron density is always recalculated from the evolving primal mass; the
commuting divergence map makes this exactly the same dual continuity update.
No second density is evolved. Both trial species partitions pass before
state, U, time, counts or diagnostics are committed. Invalid trials raise.

The existing total-energy flux and momentum use pi+pe exactly once. Equal
gammas keep the ideal total pressure independent of the partition. Compression
enters pe through changing n, rather than an added source or energy reset.
The ideal common-velocity total flux is the standard MHD flux, including both
species enthalpies. There is no separate electron-current enthalpy flux here
because this bounded reduction sets ve=vi. Numerical thermal partition error
is explicit, with all residual internal energy recovered as ion pressure.

This path requires eta=0 until physical heating partition is implemented.
No Hall, electron pressure battery, collisions, conduction, electron inertia,
vacuum or formation boundary is implied. Slip impermeable conducting radial
wall, periodic z, positive fluid in every cell. It does not combine the
incompatible boundaries of `ElectronPressure` and `HallMHD`.

`species_rhs` returns the fluid derivative, Udot, heat and limiter statistics;
inherited `rhs` remains the total ideal fluid derivative. `step` checks the
fluid CFL and species positivity at both stages; it does not promise positivity
for every state satisfying that CFL or perform a hidden retry. `advance` and
ordinary fluid diagnostics are inherited. Additional diagnostics distinguish
conserved proxy from physical logarithmic entropy and separate species energy.

Independent contact and nonlinear simple-wave evolution verify moving thermal
transport/compression. They do not establish shock partition or warm Hall
accuracy; see `WARM_ADIABATIC_NUMERICS.md`.
