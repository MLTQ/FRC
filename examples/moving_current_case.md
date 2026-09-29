# Moving-current manufactured case initialization and norms

Initializes nodal positive species data and regular bulk velocity from the
independent moving-species reference, exact corner psi, and independently
integrated annular g means. n_t and explicit electric forcing come from the
analytic density derivative and parent species forces. The forcing is not fitted
to a discrete residual, so spatial truncation errors remain observable.

Physical norms use the independently coded r,z interpolation in
current_thermal_norms, not production sample maps. Compare all components of
Bdot, Jdot and full E with analytic time derivatives/fields, keeping volume,
first-annulus and maximum errors. Independent quadrature orders measure norm
uncertainty. No time integration or self-consistent bulk motion is claimed.
