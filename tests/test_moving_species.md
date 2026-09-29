# Moving parent-species verification

Uses independent analytic nonuniform axisymmetric fields with moving density,
rotation/shear and both current components. Checks exact species continuity,
mass-weighted bulk/drift energy and full flux against direct species sums,
and the reduced momentum/Ohm residual against independently differentiated
parent species forces. Both density parts of w_t are retained explicitly.

The local energy law is checked with fourth-order physical-space differentiation
of direct species flux and observed refinement of the residual. Faraday is
checked independently against analytic time derivatives. Manufactured force and
heat are exposed, so none of these is claimed to be an unforced formation run.
