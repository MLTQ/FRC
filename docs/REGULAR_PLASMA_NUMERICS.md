# Bounded moving two-species assembly

This joins the verified regular fluid/current/species components without changing
their archived implementations. It advances the actual regular point state
`(n,n vr/r,n vtheta/r,n vz,n se,n si)` and harmonic magnetic increments together.
It is a filled, axisymmetric, periodic cylinder with impermeable radial numerical
fluid flux, fixed poloidal flux at axis/wall and the existing compatible natural
toroidal boundary. It does not represent the proposed quartz tube's vacuum,
ends, material response, coil field, ionization or transport losses.

## Original equations and energy

Let M=mi+me, mu=mi me/M, alpha=me/M, beta=mi/M, delta=beta-alpha,
w=J/(e n), ve=v-beta w, vi=v+alpha w. Physical electron and ion masses are kept.
The unforced collisionless quasineutral parent equations are

```
n_t + div(n v) = 0
M n (v_t + v.grad(v)) = -grad(pe+pi) + J cross B - div(mu n w w)
(n ss)_t + div(n ss v) = -as J.grad(ss)/e;  ae=-beta, ai=alpha
E = -v cross B + delta w cross B + (alpha grad(pi)-beta grad(pe))/(e n)
    + (mu/e)[w_t + v.grad(w) + w.grad(v) - delta w.grad(w)]
B_t = -curl(E),  mu0 J = curl(B),  div(B)=0
ps=(gamma_s-1)n hs, hs=e exp(ss)(n/nref)^(gamma_s-1)/(gamma_s-1).
```

Their physical energy is M n |v|²/2 + n(he+hi) + mu n |w|²/2 + |B|²/(2mu0).
The original local energy flux is the sum of each species' kinetic-plus-enthalpy
flux and E cross B/mu0. The current implementation's global work identities do
not prove this local flux across every reconstructed element interface.

## Discrete assembly

RegularCoordinateFluid supplies one actual `n_t`. RegularCurrent uses that same
density and derivative for physical current inertia, ideal/Hall and Lorentz work.
CurrentTransport adds current advection/stretching/self drift and its paired
bulk drift stress. RegularSpecies adds the two entropy-density updates and
pressure forcing. No independent density update, temperature reset, atmosphere,
global energy subtraction, limiter or added heat is present.

The species potential is analytically in the compatible curl kernel: only its
reduced forcing enters induction. Its gradient is restored to full E for the
original Ohm equation. KernelInduction's strict small-signal residual checks
are unchanged. A localized near-axis toroidal magnetic perturbation passes its
inner physical norm but fails its outer coefficient backward norm. Captured input
disproved the initial pressure-cancellation hypothesis. That failure must remain
visible, not bypassed by loosening tolerance. Geometry caching shares only fixed
maps/factors and is bitwise checked
against fresh operators after changing density and both entropies.

The assembled semidiscrete energy uses the actual grid fluid EOS plus physical
quadrature current and magnetic energy. A separate root-polynomial/Fourier audit
reconstructs the continuous physical energy and its derivative. These measures
are reported separately. Shared centered fluxes introduce no explicit numerical
diffusion; interpolation, dispersion and RK4 truncation remain finite errors.

## Time and spectral verification

Classical RK4 advances all state components and magnetic increments in the same
four stages. Each stage validates the physical state; returned point and sampled
densities must be positive. Failure is atomic. There is no automatic CFL claim,
adaptive retry or energy correction. In particular fixed-current entropy tests
do not determine the time step for the full moving current system.

The continuum wave is linear, with physical inertia, distinct species gammas,
radial Bessel structure and axial Fourier phase. Initial conditions and analytic
rates are independent; no discrete eigenmode or source is fitted. Nonlinear
amplitude, temporal, mesh and physical quadrature controls remain explicit.
Contact has variable n, complementary pressures, common axial speed and zero J;
its full E is nonzero and curl-free. Instantaneous contact agreement is not a
completed transit test.

A separate full stationary linearization uses n_t=-div(n0 v), the physical
pressure perturbations, Lorentz and ideal/Hall/inertial induction. Nonlinear
current transport and drift stress are quadratic at rest. The actual nonlinear
driver must converge to this operator under mixed perturbations. All active
degrees and branches are retained; the full physical energy Hessian is the
positive perturbation metric. An imaginary stationary spectrum supports only
the measured frozen problem, not arbitrary nonlinear or three-dimensional
plasma stability or particle confinement.
