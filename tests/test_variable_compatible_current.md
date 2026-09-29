# Variable-density compatible coupling verification

Tests actual annular density polynomial means and adjoint work, immutable input
and rejection of a positive-cell/negative-reconstruction state; SPD inertia,
rejection of reconstruction overflow and unrepresentable inertia, physical
current energy; second-order central directional differences of the
local changing-density/current energy; retained local harmonic cross work.

An independent direct solve of the original three-block variable-density
equations checks the unchanged induction algorithm. A physical polynomial
manufactured E, analytic curl-B current and B_t verifies the moving-density
forcing without using the production operator to construct its source. The
global current/magnetic exchange includes the density term explicitly.
The nonpolynomial continuum study reference also has an independent complex-step
check of curl B, Faraday and the prescribed particle-flux continuity identity.

Charge fluxes must pair with the physical H/M and cancel graph divergence.
Both species gammas verify lumped pressure work, entropy and U conservation
separately. This graph check is not a physical high-order pressure test; formal
continuum convergence is recorded separately. No formation or full time-step
claim follows from these instantaneous bounded component tests.
