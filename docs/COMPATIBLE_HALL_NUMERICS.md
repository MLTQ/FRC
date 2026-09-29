# Compatible compressible Hall reference

This is a bounded, fully filled, periodic cylindrical reference with a conducting
radial slip wall. `HallMHD` advances mass, all three momentum components, total
energy and all three magnetic components. It does not replace the formation
interface or supply a plasma-vacuum transition. The original `ToroidalMHD`
default and failed experimental `HallCurl` remain available for reproducibility.

## Declared continuum reduction

For fully ionized singly charged ions with prescribed mass mi, ne=rho/mi,
massless electrons and pe=0:

    ∂t rho + div(rho v) = 0
    ∂t(rho v) + div[rho vv + (p+B²/(2mu0)) I - BB/mu0] = 0
    E = -v cross B + eta J + (J cross B)/(ne e)
    J = curl B/mu0,  ∂t B = -curl E,  div B = 0
    ∂t U + div[(K+epsilon+p) v + E cross B/mu0] = 0
    U = K+epsilon+B²/(2mu0), epsilon=p/(gamma-1).

The full Maxwell stress enters momentum once. Hall adds neither a second force
nor Joule heating: J dot (J cross B)=0. Prescribed eta produces eta J²; it is
diagnosed once, while total energy receives only flux divergence. The current
reduced code assigns the resulting internal energy to the represented gas.
It does **not** solve electron/ion heating partition or justify cold physical
electrons under finite resistive heating. Independent separate thermal equations,
transport/ionization and validity checks remain required for the target regime.

The axis is regular, psi=0; the conducting wall has constant psi and zero
tangential electric field. Axial endpoints are periodic. Positive carriers are
an interior assumption, not a vacuum floor. A failed physical or nonlinear
trial raises without changing accepted state, time or physical accounting.

## Why the earlier axis operator failed

The older map differentiated an energy based on averaged poloidal fields, while
its Etheta interpolation was linear in r. Algebraic adjoint work cancellation
did not make its physical electric curl equal its magnetic coefficient update.
The first-annulus gdot error approached 21.46%, even while volume errors fell.
HALL_WORK_VERIFICATION retains that failure and its original source fingerprint.

Use x=r²/2 with volume element 2pi dx dz. Interpolate psi bilinearly in x,z:

    Br = -psi_z/r,  Bz = psi_x,  Btheta = r g.

Choose electric fields in the corresponding spaces:

    Er = r f,       f: constant in x, linear in z
    Etheta = chi/r, chi: linear in x and z, zero at axis and wall
    Ez = h,        h: linear in x, constant in z.

Then Faraday is exactly

    psi_dot = -chi,    g_dot = h_x - f_z.

Differentiating these physical fields agrees with the coefficient Faraday map C
inside every element, including regular Br,Etheta proportional to r at the
axis. There is no fitted first-cell multiplier. The compatible poloidal energy
is the positive volume integral of this field, explicitly selected by the
state's `magnetic_quadrature='compatible'` metadata. Existing states are not
silently converted or adjusted to recover their old pressure.

The radial Br basis on [a,b] is a(1-xi)/r and b xi/r, where
xi=(r²-a²)/(b²-a²). Its exact annular mass moments are integrated with fixed
16-point radial Gauss quadrature; the axial quadratic integral is exact.
This quadrature is independent of the common electric/Hall quadrature (default
5 points per direction), and audited with independent 24-point physical energy
integration. Toroidal magnetic/rotation energy uses exact mean r²=(a²+b²)/2.

## Current, skew Hall work and local transfer

Let x now denote the assembled magnetic coefficient vector (psi,g), H its
symmetric energy Hessian, M the positive electric mass and I its physical
sampling map. The weak work current is

    M j = -C.T H x.

Common physical quadrature projects eta J and J cross B/(ne e). The resistance
matrix R=I.T W eta I is positive semidefinite. The Hall matrix X, including
variable ne and the actual midpoint magnetic field, is exactly skew because
its off-diagonal component blocks are negative transposes. Thus

    M e = (X+R) j,
    dHmag/dt = -j.T (X+R) j = -j.T R j.

This is a global assembled identity. Projected local J dot E need not vanish
pointwise. Face Poynting quadrature separately approximates E cross B/mu0;
its radial/axial trace products use the same compatible axis moments. Its
local difference from magnetic energy rate is measured, not inserted as heat
or used to repair energy. That difference converges at about second order in
the nonlinear snapshot. The physical field/current norms remain partly first
order because the lowest-order magnetic spaces are piecewise in one direction.

The construction follows the compatible weak-curl/skew-work principle discussed
in [Laakmann et al., structure-preserving Hall MHD](https://arxiv.org/pdf/2202.11586).
That incompressible analysis does not verify this cylindrical compressible
implementation; independent polynomial, continuum-wave and nonlinear tests do.

## Nonlinear midpoint and fluid integration

At fixed carrier density within the magnetic split substep, solve

    M j + (dt/2) C.T H C e = -C.T H x_old
    M e = (X(B_mid)+R) j
    x_new = x_old + dt C e,   B_mid = B((x_old+x_new)/2).

Picard iteration updates B_mid. Diagonal mass and physical impedance scaling
balance the mixed system. Cached sparse LU is only a preconditioner for the
actual new linear system; bounded GMRES checks its residual and falls back to
a fresh factorization when needed. Both nonlinear Ohm and linear residuals
must pass. A rejected trial may alter the computational cache, not the accepted
physical state. No artificial resistivity or field damping is added to converge.
Quadratic magnetic energy changes by -dt j.T R j at convergence; the positive
heat diagnostic uses the same midpoint J. The total fluid energy update is
the midpoint Poynting divergence and has a separate local truncation error.

A full step is ideal-fluid SSPRK2(dt/2), Hall/resistive midpoint(dt), then
ideal-fluid SSPRK2(dt/2). The positive limited fluid reconstruction supplies
actual carrier density at common quadrature points. Momentum is unchanged only
inside the magnetic substep, and evolves in both ideal substeps. The original
Rusanov flux, UCT and reconstruction limiting retain their declared numerical
diffusion. Physical Ohmic heat excludes that numerical dissipation. There is
no complete-update positivity theorem or shock accuracy claim for this split.

The time limit controls explicit fluid wave stability. Implicit midpoint avoids
an explicit whistler stability restriction, but does not bound its phase error.
Users must supply time caps and demonstrate refinement. Nonlinear convergence
or the post-Hall pressure/fluid step can require a smaller step; automatic
subdivision is not implemented. Failed steps are atomic. `hall=False` retains
resistance; with eta=0 its RHS is the compatible ideal solver and its complete
step equals two ideal half steps.

## Evidence and next gate

COMPATIBLE_HALL_VERIFICATION contains actual coupled wave propagation, signed
propagation/amplitude/time/tolerance controls, finite-amplitude variable-density
evolution and independent resistive field/heat decay. The legacy tests remain
passing. No evidence here establishes formation, vacuum coupling, electron
inertia, finite electron pressure, real ends, three-dimensional stability or
particle confinement. TWO_TEMPERATURE_PLAN defines the next bounded closure.
