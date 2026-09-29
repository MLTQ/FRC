# Higher-order compatible fields and exact harmonic storage

This component addresses the diagnosed slow-mode cancellation between bulk and
Hall electric terms. The physical finite-mass equations and original parameters
are unchanged. It is currently a constant-background **linear verification
component**, not a replacement for the nonlinear moving finite-mass solver or
an accurate FRC formation simulation.

## Spatial field definition

Use x=r²/2. On every physical interval, a local degree-p polynomial interpolates
the nodal values, using p+1 interior/periodic nodes including both endpoints.
Degrees p=1,3,5 are explicit. Interior stencils overlap; at the radial ends they
remain within the physical domain. This is our specific local construction,
not an assertion that a published spectral-element convergence theorem applies.

Let N_i be the nodal functions, h_j the interval lengths and

    D_j f = (f[j+1]−f[j])/h_j,
    E_j = −h_j sum_{i<=j} dN_i/dx.

Then d(N f)/dx=E D f, and the integral of E_j over interval k is h_k delta_jk.
Edge coefficients therefore represent actual interval averages. The analogous
construction in periodic z uses unwrapped local nodes with wrapped coefficient
indices. Derivatives/edge functions may jump between physical cells; nodal
functions remain continuous. No nonlinear limiter or positivity property follows.

Nodal/edge differential-form interpolation and commuting projections follow
the general framework discussed by [Kreeft, Palha and Gerritsma](https://arxiv.org/abs/1111.4304).
Compatibility alone is insufficient to prove stability: a bounded cochain
projection is an additional hypothesis in the analysis of
[Arnold, Falk and Winther](https://arxiv.org/abs/0906.4325). We measure our
particular overlapping basis rather than claim that theorem without proof.

With tensor products understood, the physical fields are

    Br = −Nx (dNz/dz) psi / r,
    Bz = (dNx/dx) Nz psi,
    Btheta = r Ex Ez g,
    Er = r Ex Nz er,
    Ez = Nx Ez ez,
    Etheta = Nx Nz (r_i et_i) / r.

The familiar incidence equations remain exactly the same mathematical map:

    psi_t = −r_i et_i,
    g_t = D_x ez − D_z er.

Thus differentiated reconstructed B equals the physical −curl(E) representation,
and C grad=0. Boundary psi rates vanish at the axis and wall. The full scalar
gradient has er=D_x phi, ez=D_z phi and et=0. Uniform periodic axial electric
voltage remains representable. The origin is regular: psi(axis)=0; an inadmissible
N_axis/r basis is excluded rather than integrated through a singularity.

## Physical work, current and inertia

Common positive physical 2 pi r dr dz Gauss weights define electric mass M and
magnetic energy Hessian H. Tensor-product Gram matrices avoid storing enormous
multi-dimensional sampling matrices. Their factors are geometry-only. The
declared quadrature order remains explicit; independent order changes and
root-product polynomial evaluation check integration and interpolation errors.

    M J = −C.T H B,
    W_B = B.T H B/2,
    I = mu M/(n e²),   mu=mi me/(mi+me),
    M E = I J_t + f,
    B_t = C E.

Original finite-inertia elimination, residual tolerances and linear refinement
are retained. The current is obtained from actual magnetic coefficients and
physical masses, never fitted to the desired wave or velocity. A separately
assembled original three-block mixed system tests the elimination. Independent
polynomial B and curl(B)/mu0 verify physical current, including wall/axis traces.

For constant n and B0, the linear forcing is the projection of
(−v+delta J/(ne)) cross B0, with delta=(mi−me)/(mi+me). The exact complete
pressure gradient is represented by a nodal potential
(alpha pi−beta pe)/(ne), alpha=me/(mi+me), beta=mi/(mi+me), and restored after
the curl solve. Full original Ohm/Faraday backward errors and rate-relative
Faraday errors are reported separately. These are linearized equations;
nonlinear current-inertia vorticity, n_t and energy/thermal terms are not omitted
from a purported nonlinear evolution—they have not yet been coupled here.

Velocity comes from independently initialized conservative cell averages.
Degree-3/5 annular polynomials r^s P(x), s=0 for rho/mz, 1 for mr and 2 for L,
are combined with periodic axial polynomials. Reconstruct rho before dividing
momentum. This does not use analytic current or substitute a reference velocity.
Nonpositive reconstructed rho rejects; no atmosphere, floor or limiter is added.

## Harmonic background and small-signal precision

High-order spatial cancellation exposes loss of digits when a tiny field is
stored or differentiated beside a large uniform B0. The initial rounded-total
study is preserved separately. Its finer-grid error is not declared spatial
convergence merely because all linear solves have small residuals.

Represent B=B0 ez+delta B with an explicit constant bias and independently
stored delta psi,g. Both radial delta-psi boundary traces vanish. This is the
same physical field in an affine coordinate representation. The constant
harmonic has exactly zero current, and

    integral B0 delta Bz dV
      = 2 pi B0 integral[delta psi(R,z)−delta psi(0,z)] dz = 0.

Ampere therefore acts on delta B, while full field samples retain B0. Global
energy keeps harmonic and increment components separately; local B² includes
the spatially nonzero cross term. No threshold removes weak modes. Tests cover
currents down to 1e-12 of the polynomial reference beside an unchanged bias,
with error measured relative to the weak current itself.

Already rounded input cannot recover lost digits. `rounded_total` is an explicit
lossy interchange format and is never silently used as an exact increment.
This harmonic decomposition follows the fixed-flux periodic-cylinder boundary;
it does not supply general equilibrium balancing, moving boundaries or spatial
coil fields. Exact wall values of the analytic Bessel mode are imposed from its
specified root boundary, not obtained by fitting velocity to numerical current.

## Required moving-coupling contract

Do not drop this geometry into the old finite-mass state as a one-line switch.
That would change physical magnetic/current energy, pressure recovery and field
traces without deriving their common work. Before adoption:

1. Define variable positive density/inertia quadrature, physical local B/J energy
   and its exact n/B derivatives using the same declared basis. Preserve the
   actual conservative mass flux and quasineutral n_t; do not evolve a second
   independently chosen density.
2. Derive consistent velocity and cylindrical vorticity, stress, Poynting and
   current/enthalpy work with those fields. Exact global cross cancellation does
   not permit omitting local harmonic cross energy or bias forces.
3. Derive electron/ion transport and the pressure gradient as an adjoint pair
   with the new current. Weak C grad=0 is not a proof that old dual-cell particle
   fluxes still match the new physical current representation.
4. Repeat original full species/manufactured equations, all wave pressure and
   magnetic rates, finite-time contacts/waves, local energy/entropy, axis/max,
   time/quadrature/mesh and atomic rejection tests. Preserve original physical
   inputs and previous failures; no fitted initial states or energy repairs.

Then reconnect the original heated plasma edge and exact exterior vacuum,
physical collisions/transport/ionization and justified boundaries before any
formation study. No physical diffusion is added by this component. It has no
shock/entropy theorem and establishes neither nonlinear stability nor confinement.
