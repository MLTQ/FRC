# Conditional Coulomb transport and local scales

`coulomb_coefficients` uses NRL Maxwellian thermal collision times in SI,
electron-ion Coulomb-log branches and same-species ion-ion logarithms. Inputs
are ne [m^-3], Te/Ti [eV], mi [kg], integer Z. Electron thermal speed exceeds
ion thermal speed in the implemented electron-ion branch. Unsupported branches
and nonpositive logarithms fail explicitly; no logarithm floor is applied.

Momentum resistivity is me nu_ei/(ne e²), the perpendicular friction coefficient
in the strongly magnetized classical limit. For Z=1 the separately named
parallel/unmagnetized estimate is smaller by 1.96. Neither is the Pedersen
conductivity, which involves Hall coupling. Nonunit-Z parallel transport is
not approximated by the same factor. Thermal mean free paths use sqrt(kT/m),
not the alternative sqrt(2 kT/m) speed convention.

`regime` evaluates positive plasma states at an explicitly supplied physical
length/time scale, field and current. It reports Hall/electron-inertia/Debye
lengths, carrier drift, Coulomb collisionality and scope-specific review flags.
It returns classical coefficients as null when the closure is unavailable.
Zero-density vacuum is excluded rather than filled; it has no plasma transport.

Thresholds 0.1, lnLambda=5 and light-transit ratio 0.01 are conservative review
criteria, not necessary/sufficient validity theorems. No green validity flag is
inferred from their absence. Neutral contributions, distributions, ionization,
thermal transport and instability-enhanced collisions remain unmodeled. Scans
do not choose laboratory parameters or automatically change prescribed eta.
Primary equations and limitations are in TRANSPORT_REGIME.md.
