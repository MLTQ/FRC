# Regular species instantaneous study and independent audit

`run_case` applies both physical species pressure/entropy terms to nonuniform
density and distinct temperatures. It records original analytic nodal rate
errors, entropy/quadratic/work budgets, and the independently solved magnetic
plus current-inertia energy response. Density, bulk momenta and angular momentum
are unchanged by this isolated pair. Point-transfer unskew/skew controls retain
the rejected handoff constructions. No tolerance is selected from their outcome.
Induction solves the analytically non-gradient forcing; its exact compatible
potential is restored to E and the original full Ohm residual is checked.

`audit` reconstructs saved rates and electric/current fields with independent
root-product polynomials and FFT coefficients at a caller-selected quadrature.
It reports whole-domain, first/last-fluid-element and maximum errors against
the original analytic pressure and entropy equations. It measures actual
continuous thermal work separately from nodal/grid work and retains strong
div(J), which need not be zero despite graph charge continuity.

These are manufactured instantaneous rates. They do not establish nonlinear
moving, shock, free-edge, physical closure, formation or confinement validity.
