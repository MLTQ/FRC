# Independent joint-transport reference and physical norms

Reuses the original two-temperature translating contact and finite-amplitude
simple wave without changing density, pressures, gammas, masses or duration.
A separate instantaneous axial Euler reference has nonuniform density, three
velocity components and independently varying species pressures. Complex-step
derivatives of those continuum fields give original conservative n/m/C and total
energy rates, without calling the numerical flux or its energy gradient.

Continuous physical norms use a declared real trigonometric interpolant of the
conserved point state, sampled independently by Fourier resampling. Physical
temperatures, pressures and kinetic/internal energy are then recomputed from
that interpolant. Denser 8x/16x integration checks quadrature convergence.
Nodal-energy conservation is never relabeled as exact continuous-energy
conservation: their difference and its change remain separately recorded.
Negative interpolated density rejects instead of using a numerical floor.
