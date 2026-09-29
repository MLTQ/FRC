# Next bounded integration: compressible toroidal degrees of freedom

The independent `LinearHallCylinder` now verifies axis regularity, Hall/finite-
mass waves and the linear work pairing. Its uniform-density incompressible
pressure projection is not a compressible solver. The next concrete task is to
extend the bounded conservative MHD core with toroidal field and rotation while
preserving the existing poloidal API and its verification records.

## Conservative variables and contracts

For ideal/massless-electron MHD, use axial angular momentum density
L=r rho vtheta and toroidal induction variable g=Btheta/r, alongside rho,
rho vr, rho vz, total energy and corner psi. At the axis both Btheta and vtheta
are regular O(r); L=O(r²) and g is finite. Conversions must integrate/reconstruct
these geometric weights correctly, not divide a cell-average L by an arbitrary
radius and claim exact kinetic energy. State layout changes need an explicit
separate type or migration; existing four-variable states must not be silently
reinterpreted.

For the poloidal vector vp and Bp, the continuum fluxes are

    d_t L + div_p [r(rho vtheta vp - Btheta Bp/mu0)] = 0
    d_t g + div_p [(Btheta vp-vtheta Bp)/r] = 0,

where div_p F=(1/r)d_r(r Fr)+d_z Fz. The second equation follows directly from
Faraday's theta component; using the ordinary cylindrical divergence on Btheta
itself would add a spurious geometric term. Radial momentum has source

    [rho vtheta² + p + B²/(2 mu0) - Btheta²/mu0]/r.

Energy includes all three velocity/field components and uses
(E+p+B²/2mu0) vp - (v dot B) Bp/mu0. Poloidal flux continues to use the shared
Etheta and exact CT divergence. Both radial axis and outer reference wall must
have explicit torque and electromagnetic work; this does not supply vacuum.

## Independent checks before Hall integration

- Rigid rotation vtheta=Omega r, constant rho, Bz and
  p=p0+rho Omega² r²/2 verifies geometric centrifugal support.
- The combined regular pinch/rotation Btheta=c r has
  p=p0+(rho Omega²/2-c²/mu0)r². Independent Maxwell derivatives give
  Jz=2c/mu0, and radial force cancels the pressure/centrifugal terms.
- Linear torsional/poloidal waves must recover the established independent
  cylindrical references in the applicable limit. Measure actual field/velocity
  errors, torque and energy; keeping divergence small is insufficient.
- Finite-amplitude advection/stress tests and radial/axial mesh/time refinement
  must expose numerical diffusion and energy drift. Preserve old benchmark
  behavior at exactly zero toroidal variables.

Only then add the Hall electromotive force through a work-paired current/curl
operator, using the new 49-cylinder / 16-parallel records as references. Variable
density makes the pressure term and current-inertia transport nontrivial. In
particular, summing parent species momenta includes a drift momentum flux
beta J J (beta=mi me/[(mi+me)n e²]); ignoring it while counting drift energy is
not a valid nonlinear finite-mass model. Thermal/current advection, geometry
and collision exchange must come from the same species reduction.

Full heated edge/vacuum coupling, spatial sources, actual ends, kinetic ordering,
shock/remapping robustness and converged formation remain in EXTENDED_PLASMA_PLAN.
The unchanged original argon target remains the objective. Bounded toroidal
verification is an integration milestone, not an alternate final simulator.
