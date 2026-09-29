# Affine harmonic magnetic storage

Stores B=B0 ez+delta B without first rounding the small perturbation into a
large total psi. B0 is an explicit constant physical axial bias; delta psi
vanishes at the axis/wall and is periodic. In this mixed bounded domain its
exact current is zero and its cross-energy with delta B integrates to zero:
integral B0 delta Bz dV=2 pi B0 integral[delta psi(R)−delta psi(0)] dz.
Consequently Ampere acts on increments, while all physical field samples and
local energy retain the full bias. Local cross-energy is not zero and is kept.

`packed_increment` is the magnetic rate/current coefficient space; `samples`
returns the complete physical B; `current` applies the actual compatible weak
Ampere operator to the increment. `energy_components` keeps harmonic, increment
and exact global cross terms explicit. `energy_density` integrates the full
local B², including cross terms. `rounded_total` is a deliberately lossy export,
not an evolution/current input. Immutable copied arrays prevent accidental
changes to the accepted representation.

There is no wave fitting, amplitude threshold, changed mass, omitted bias force,
field clamp or energy reset. A rounded total input cannot recover digits already
lost; callers must initialize/store increments directly. This exact constant
harmonic is specific to fixed-flux periodic-cylinder boundaries, not a general
equilibrium decomposition, arbitrary coil field or moving vacuum interface.
Nonlinear local work/energy/pressure and evolving boundary coupling still need
derivation before using this representation in the formation integrator.
