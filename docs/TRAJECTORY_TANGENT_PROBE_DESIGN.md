# Next spatial probe after the complete rejection replay

Status: design only, not implemented or executed. First collect the live
rejection campaign and verify its entire accepted-prefix/failure protocol.
Do not run this large calculation concurrently with that worker.

The density observations locate a developed low-axial-mode departure spread
across radius, with internal trace separation. They do not identify a spatial
eigenmode, establish an eigenvalue or distinguish spatial versus temporal error.
The next bounded diagnostic should directly apply the unchanged full spatial
RHS to a direction derived from the actual early trajectory. A dense22,208 by
22,208 Jacobian is unnecessary for this question.

Use the exact matched fixture F8/64,B32/64 and actual split inverse. Preserve
the original physical masses, density, both entropy equations, current stress,
wall space, bias, tolerances and induction residuals. Explicitly retain the
separate species-pressure gradients; the total pressure alone is constant.

Let U0 and dU be the original conservative initialization and independently
defined linear tangent, and let b0=0 and db be the original initialized
admissible magnetic tangent. Using the saved step10 (16.667us), define

    eU = [U10 - Uinitial] - epsilon*expm1(lambda*t10)*dU
    eb = [b10 - binitial] - epsilon*expm1(lambda*t10)*db.

Keep the entire fluid/current/magnetic error direction. Do not retain only
selected axial modes, radial elements or density. Hash-check the original
input, every selected state and the independent continuum reference. Record
the exact initialized tangent, native error coefficients and their coordinate
representation, including pack/unpack differences and constraint residuals.
The direction is a diagnostic input; no evolved state is replaced by it.

Use ModalWallCoordinates' fixed physical scales and admissible tangent space.
Normalize the packed direction by its infinity norm, recording that factor.
Measure the actual centered derivative of the full original RHS about U0:

    D_h F(U0)[v] = [F(U0 + h*v) - F(U0 - h*v)] / (2*h).

Use h=1e-4,1e-5,1e-6 as separately recorded amplitude controls; preserve any
rejection. Apply the same protocol to the original independent continuum
tangent as a control. Save both signs' full fluid/magnetic rates, solved E,
current diagnostics, positivity and all original residuals, plus the base rate.
The nonzero discrete base-stationarity defect must be measured explicitly.

Coordinate projection of D_h F onto v and its orthogonal remainder can be
reported with the declared fixed metric. Such a Rayleigh quotient is only a
directional diagnostic, not an eigenvalue unless the full residual supports
that claim. Independently reconstruct physical field directions and their
rates at quadrature16/20; retain separate component norms, both pressures and
current/E rather than using only a mixed-coordinate norm. The known-mode
control must retain its comparison against the independent continuum rates.

Before causal interpretation, require amplitude agreement, original strict
equation checks, full fresh replay of all raw data/metrics and physical
reconstruction checks. A spatially growing direction still requires continuum
and numerical comparisons to determine whether its growth is physical. This
probe alone does not replace full opposite-sign/time/mesh controls or establish
the cause of the trajectory failure. No fix is selected in advance; no filter,
diffusion, floor, artificial mass, pressure/energy repair or relaxed bound is
authorized by this diagnostic. Formation and the full sandbox gates stay open.
