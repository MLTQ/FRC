# Frozen-material electron current inertia

`MaterialCurrent` advances the linear toroidal current response with physical
electron inertia and scalar momentum resistivity. It reuses the verified
full-space magnetic operator and plasma-only current mass. Vacuum remains an
elliptic magnetic region with exactly no carrier mass, heating or resistance.
It is separate from the moving solver: no incomplete moving electron-energy
or Hall model is silently introduced.

At physical quadrature points, density must be finite and positive inside;
resistivity is finite and nonnegative. Scalar, per-cell or callable profiles
are sampled once and frozen. `electron_density=None` explicitly selects the
zero-inertia model comparison. No density floor or anomalous resistivity exists.

The mixed unknowns are induced a and w=mu0 Jtheta/r. Ampere's law and midpoint
Ohm's law form a sparse block solve. Positive weighted masses integrate de²
and eta/mu0. It counts full-space induced magnetic energy and electron drift
kinetic energy separately, plus positive per-cell collisional heat and source
work. No residual is redistributed. Its derivation is CURRENT_INERTIA.md.

`step` is pure with respect to its input field and accepted state. Returned
`CurrentUpdate` includes both endpoint energies, current, heat and separate
Ampere/Ohm/vacuum residuals. The block solve advances increments, preserving
stationary current exactly and avoiding subtraction of large endpoint states.
Ohm residual is the backward error of that increment equation, including its
right-hand side. Diagonal congruence scaling conditions the mixed linear solve at small axis
weights; it changes neither equations nor physical coefficients. A mass/length
scale handles exactly zero inertia/resistance blocks without modifying density.
The exactly massless, zero-resistance limit delegates to the existing exact
frozen-flux branch, avoiding cancellation of a tiny material coefficient when
the prescribed bias is much larger. Its heat and current kinetic energy are zero.
The Ampere residual is still explicitly evaluated in this exact-flux branch.
Implicit midpoint is second order but not L-stable; resolved time integration
remains required. Coefficients/geometry are immutable
for one operator. Rebuild after any change; do not mutate its stored matrices.

The component omits material motion, electron pressure, Hall feedback, heat
transport, ionization and displacement current. It is a controlled linear
induction limit, not a full plasma-vacuum edge solution or a formation solver.
