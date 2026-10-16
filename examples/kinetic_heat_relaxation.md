# Complete linear kinetic temperature-wave relaxation

Runs `python -m examples.kinetic_heat_relaxation OUTPUT` using the independent
`kinetic_heat` component. Uniform, stationary, infinitely massive ions with Z=1;
small electron temperature Fourier perturbation, zero density/current response,
constant positive Coulomb logarithm absorbed into nu0. Space is periodic. There
is no plasma-vacuum interface or lab formation prediction in this fixture.

Degree8/12/16 velocity expansions evolve full temperature/heat histories at
K=k*vT/nu0 =0.001,0.03,0.3,1. Full electron-electron linearized Landau terms
and velocity-dependent electron-ion angular scattering remain enabled. Matrix
exponentials eliminate a time-discretization parameter; velocity truncation and
the eigenvector condition remain measured. No state clipping or fitted damping.
Complete coefficients, operators, histories, physical moment diagnostics and
source identity are saved. The example adds no unit tests or test framework.

Local comparisons: solve the independent steady, zero-current Lorentz equation
analytically in speed. Its conductivity is <s^7>/6=64*sqrt(2/pi), using vT=sqrt(T/m).
For the full Z=1 collision operator compare the long-wavelength history with
Braginskii's approximate coefficient3.16 in units n*vT²*tau_e,
tau_e*nu0=3*sqrt(pi/2). This published coefficient is approximate, not an exact
target to fit. The run also reports its own converged steady collision solve.
Source: [Kupfer, Harvey et al. (1996), classical comparison section](https://compxco.com/Kupfer_Harvey_etal_1996_PhysPlasmas_3_3644.pdf).

Collisionless reference: independently integrate Maxwellian characteristics
and enforce zero current. With u=k*vT*t, the normalized electric response A
satisfies a second-kind Volterra equation

    A(u) + integral_0^u A(w) [(u-w)^3-3(u-w)] exp(-(u-w)^2/2) dw
         = [1-5*u²/2+u^4/2] exp(-u²/2).

The temperature amplitude is

    [1-7*u²/6+u^4/6] exp(-u²/2)
      + integral_0^u A(w) [5(u-w)-(u-w)^3]/3 exp(-(u-w)^2/2) dw.

These formulas follow direct Gaussian density, momentum and energy integrals;
they do not use production moments or its streaming matrix. Trapezoidal
Volterra integration at2000/4000 intervals measures reference integration error.
Compare complete0<=u<=4 histories with collisionless degree12/20/28 expansions.
Finite-degree recurrence beyond the measured interval is not excluded.

Use `--degrees 20 24 28` to refine the velocity expansion without altering the
physical equations or time interval. `--reference-intervals 16000` refines only
the independent characteristic integration (and retains its half-resolution
comparison). These choices are recorded explicitly; prior outputs remain frozen.

The kinetic, zero-current heat-wave setup and the limitations of local/steady
closures are discussed in [Brodrick et al.](https://arxiv.org/abs/1704.08963).
The actual evolving pulse has finite drift, nonuniform B, ion dynamics and an
edge; adoption requires those physics and boundary budgets. This baseline does
not justify replacing them with the uniform linear calculation.
