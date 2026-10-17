# Conservative velocity shape and explicit mixing heat

The explicit research option `projection_method='variance'` addresses local
cooling in U=E-K and density-weighted velocity extrapolation. All occupied source
points, fine cells, old velocity knots, physical masses and populations remain.
It is not adopted formation physics and uses no floor or cutoff.

Let rho_s,u,U_s denote the quadrature source, rho_t its cell-average target
density, M the common block mass and c=P/M its physical mean velocity. Fit
source velocity in occupied-volume L2 with the target mass mean constrained to c.
Write the zero-mean fitted shape as w. Transfer

    v = c + alpha w,
    U_t,c = integral_c [U_s + rho_s |u-v|²/2] dx,
    alpha = 2 integral rho_s w·(u-c) dx /
                  [integral rho_s |w|² dx + integral rho_t |w|² dx].

Both alpha=0 and the displayed nonzero root conserve energy. Choose whichever
has smaller occupied-volume velocity error; ties choose the constant mean.
If w is exactly zero, v=c. No near-zero threshold or clipping of alpha is used.
The factor is solved as part of the velocity/moment projection before assigning
heat; it is not a correction of a measured global energy residual. Source and
target density forms both matter because source density can vary within a cell.
Target kinetic energy uses exact P1 integrals; source moments retain the original
Gauss points and weights.

Because integral rho_t w=0 and P=M c, the summed energy difference is

    delta E = alpha integral rho_s w·(c-u)
            + alpha²/2 [integral rho_s |w|² + integral rho_t |w|²] = 0.

Momentum is conserved independently. Direct source heat plus squared velocity
error is nonnegative for admissible source heat. Source heat is still recovered
from float64 source total moments; negative values reject. Uniform velocity
shifts leave the mathematical projection invariant. Observed relative momentum
and energy residuals above 1e-12 reject without retry or energy reset.
An exactly cold stationary block has zero energy scale and exactly zero
redistribution; its residual is recorded as zero without inventing a scale.

Local total energy changes with the transferred momentum. The caller records
its redistribution L1 and mixing heat as numerical errors, not physical heating.
A single shape can still be inaccurate and its factor may be far from one or
negative. Selecting the constant root erases block velocity structure, so that
choice and both factors are reported explicitly. Measure local errors, zero-source evolution and full source/grid
convergence; positivity and conservation alone do not establish accuracy.

This derivation is specific to this code. Related staggered remapping work:
[Kenamond, Kuzmin and Shashkov (2021)](https://doi.org/10.1016/j.jcp.2021.110254).
This is not their algorithm.

Inputs are one species' flattened source quadrature fields, block-local cell
indices, model, normalized physical momentum and cell energies. Returns P1
velocity, P0 heat, signed local energy changes and diagnostics.
