# Finite-mass compressible wave oracle

Linearizes the exact two parent species about uniform positive n, pe, pi,
zero flow and axial B0. rho=(mi+me)n, cs²=(gamma_e pe+gamma_i pi)/rho,
delta=(mi-me)/(mi+me), mu=mi me/(mi+me), Lambda=1+mu K²/(mu0 n e²).
The momentum/pressure amplitudes retain the original continuum equations;
both magnetic rates acquire Lambda and Hall's coefficient acquires delta.

In equal-energy coordinates, the six-dimensional Hermitian problem is the
independent compressible generator with vA/sqrt(Lambda), delta DH/Lambda.
Magnetic amplitudes are converted back by sqrt(Lambda). The independently
eliminated dispersion relation is

    (Lambda omega²-k² vA²)
    [Lambda omega⁴-K²(Lambda cs²+vA²)omega²+k²K²vA²cs²]
      = omega² delta² DH² k²K² (omega²-K²cs²).

Full E includes both species pressure gradients and mu J_t/(ne²). Species
velocities are v+alpha J/(ne), v-beta J/(ne); perturbation pressure fractions
are gamma_s ps/(gamma_e pe+gamma_i pi), including unequal gammas. The three
branches are unforced compressible modes, not a prescribed bulk constraint.
Physical Maxwell/species differential residuals must verify the oracle.
derivatives differentiates the stated Bessel/Fourier spatial basis analytically,
independently of the amplitude equations. This avoids finite-difference loss of
digits in the slow mode's curl of a much larger nearly curl-free electric field.

The inherited amplitude is the equal-energy norm [m/s], now including current
kinetic energy. The inherited gamma is gamma_e. J1(qR)=0 ensures the bounded
slip/zero-Btheta/fixed-psi traces; full Ez is free. These are linear solutions
with quadratic finite-amplitude residuals, not exact nonlinear waves or vacuum
formation fields. Explicit me=0 and physical equal masses are limiting controls.
