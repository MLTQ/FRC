# Moving resistive plasma and physical heat

`ResistiveMovingPlasma` extends the ideal [material/vacuum reference](MOVING_PLASMA.md)
with scalar prescribed electrical resistivity eta>=0, induced-current diffusion
and positive local Ohmic heating. It retains a real moving plasma surface,
magnetic vacuum to infinity and no exterior gas. It is still a verification
solver: transport/Hall validity, spatial coil sources, shocks, remapping and
physically converged formation remain unfinished. No Stage 2 gate is removed.

The underlying magnetoquasistatic relations are Faraday, Ampere and interior
Ohm law; see [MIT, magnetic diffusion laws](https://web.mit.edu/6.013_book/www/chapter10/10.5.html).
With the repository's flux-per-radian convention, inside stationary material

```
psi_t = -eta r Jtheta,
a_t = -eta Jtheta/r,    a=psi/r²,
Jtheta = -(r/mu0)(a_rr+3a_r/r+a_zz).
```

The material derivative adds the existing ideal advection rule during motion.
Uniform external bias B(t) contributes a_external=B/2. Its self energy at
infinity is excluded from the finite response energy. This is not a circuit
model or a calculation of delivered hardware energy.

## Compatible current and positive cell heating

At fixed material coordinates assemble the consistent positive plasma mass

```
M_ij = integral_plasma r³ N_i N_j dr dz,
K_ij = integral_allspace r³ grad N_i dot grad N_j dr dz,
```

with K's exactly integrated infinite exterior modes from the ideal reference.
Vacuum nodes and modes have zero material mass. They satisfy `(K u)_vacuum=0`
at every completion, where u=a_ind. The weak current is Jtheta=r q, defined by
`M q=(K u)_plasma/mu0`. There is no exterior Ohm law, conductivity or heat array.
This current includes finite-element interface derivative mismatch; its
continuum profile and boundary accuracy are independently tested.

For constant eta, the resistive evolution is
`M a_material_dot=-(eta/mu0)(K u)_plasma`. Use implicit midpoint at fixed geometry
with exact prescribed bias increment Delta B over the substep:

```
[Mfull + dt eta K/(2mu0)] u_new
  = [Mfull - dt eta K/(2mu0)] u_old - Mfull (Delta B/2),
(K u_new)_vacuum=0.
```

Mfull has no exterior rows. The bias vector only needs material entries because
of those zero rows. An eta=0 branch completes vacuum directly without solving
a singular mass-only exterior. It copies material coefficients exactly, even
when they are much smaller than the external bias.

At the midpoint define q from the consistent mass solve. The actual cell heat
and prescribed work are

```
Q_c = 2pi eta dt integral_cell r³ q_mid² dr dz >= 0,
W_source = -(pi/mu0) Delta B sum_plasma (K u_mid)_i.
```

Symmetry, vacuum stationarity and the midpoint identity imply
`E_ind_new-E_ind_old + sum Q_c - W_source = 0` up to solve/roundoff error.
Heat is integrated from current, not assigned from the observed energy residual.
It is deposited only into material internal energy. Independently integrating
physical B², including beyond-grid energy, checks this budget without relying
solely on the quadratic stiffness expression. Profile/heat convergence against
a continuum eigenmode is separate evidence.

## Why diffusion needs an accuracy controller

For a decaying eigenmode `y'=-lambda y`, raw midpoint has amplification
`(1-lambda dt/2)/(1+lambda dt/2)`. It is energy dissipative but changes sign for
lambda dt>2 and tends to -1 for very stiff modes. Exact energy accounting does
not prevent false magnetic reversal. The raw component remains available for
explicit verification, while the moving solver uses adaptive substeps.

For each trial, compare one full midpoint step with two half steps in three
independent norms: full-space induced magnetic-energy norm, physical-current
weighted norm and volume-weighted local heat-density norm. Accept the two-half
solution only when every relative estimate meets the declared tolerance
(default 1e-4). Sum its positive physical heat and actual source work; there is
no extrapolation, energy correction or new diffusivity. Rejected trials change
no accepted plasma state. The square-root step controller is bounded and all
accepted/rejected counts and maximum accepted error estimates are recorded.

Error scales retain the largest accepted amplitude over the requested span.
The initial induced energy divided by square-root material volume supplies
the local-heat norm's absolute scale. These are error normalizations, not
density/pressure/field floors. They allow negligible decayed tails to stop
dominating work. A local tolerance is not a global error guarantee, and fields
below numerical uncertainty cannot support topology or sign claims.

Trial exhaustion or timestep underflow raises atomically. Finishing on the
last permitted trial succeeds. The generic fixed-geometry wrapper interpolates
bias linearly across its requested span; nonlinear source resolution remains
the caller's outer-time responsibility. Moving resistive stages hold bias
fixed, as described next.

## Composition with actual material motion

Advance a full moving step with resistive half-step, ideal Verlet full-step,
resistive half-step. In extended time, only the ideal stage advances the clock
and uniform external bias. Each resistive stage uses the coordinates and bias
at its own fixed stage time. This symmetric Strang composition does not drive
the source twice. In the full moving solver all external-source work is thus
accounted by the ideal stage's existing time quadrature.

For each resistive cell heat Q at fixed volume, update the adiabatic invariant
by `kappa_new=kappa_old+(gamma-1) Q V^(gamma-1)`. Invert the declared ideal
geometry/axis rule to express the new material a in reference flux labels.
Labels and entropy now change physically; the subsequent ideal stage preserves
their new values. Cell and nodal inertial masses remain fixed. This is a single
gas-energy equation, not separate electron/ion temperatures.

All stages run on a copy of mutable state, labels, entropy, heat and scalar
accounting. Commit only after the whole step succeeds. A failure in the second
resistive half-step leaves the original accepted state and every diagnostic
counter unchanged. The full response-energy residual is not forced to zero:
ideal time integration, source quadrature and splitting errors must converge.

The material wave/motion cap remains an outer-time estimate from the accepted
state, with an explicit source cap for a dynamic waveform. It is not an
unrestricted positivity or mesh-validity theorem after arbitrarily strong
heating. Geometry distortion, cold boundary accuracy and stiff thermal response
remain subjects for independent dynamic tests.

## Independent continuum and time references

The [spherical decay derivation](MOVING_RESISTIVE_PLAN.md) supplies
`a=3A j1(ks)/(ks) exp(-D k²t)` inside a fixed conducting sphere, k=n pi/R,
with a continuously matched exterior dipole. It predicts full profiles,
one-sided boundary current, local time-integrated heat and exact full-space
energy. The sphere is held fixed for this component check; its magnetic force
is not balanced by a physical fluid pressure. It is not an evolving equilibrium.

Separately diagonalizing the fixed spatial operator supplies an exact
semidiscrete time target. It isolates the raw midpoint false reversal and
adaptive tolerance behavior without confusing spatial eigenvalue error with
temporal error. Actual moving FRC/bias data then checks splitting/time refinement,
material mass, physical heat, full-space response work and source-driven mesh
sensitivity. Initialized reversed flux is not formation, and positive pressure
or a small budget error is not physical accuracy.

## Remaining physical and numerical scope

Eta is prescribed and spatially uniform. No species/temperature closure,
ionization, electron/ion exchange, neutrals, radiation, Hall term or conduction
is inferred. A variable eta would require the consistent resistance matrix
derived in the precursor plan; multiplying nodal current by an arbitrary
variable coefficient has not been substituted. Material-interface pressure and
density behavior, including heated low-pressure surfaces, still require refined
free-boundary/transport validation. The inherited Q1-gradient field accuracy,
axis flux interpolation, matching dependence and possible mesh low-energy
modes remain explicit. There is no shock viscosity, numerical resistivity,
remapping, external fluid atmosphere or thermodynamic repair.
