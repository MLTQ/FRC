# Explicit graph pressure/work bridge

`CompatibleSpeciesPressure` retains the previously derived lumped nodal thermal
definition U=n w, p=(n/nref)^(gamma−1) U, A=(n/nref)^(gamma−1)/(gamma−1),
W_s=sum V_dual A U. At fixed supplied nodal density it pairs

    U_t=div_graph(w_face I/e), E_s=−grad(p/n)/e−w_face grad(A)/e,

using the **new Hessian-derived** compatible charge I, not the old g stream.
The reciprocal logarithmic entropy mean is unchanged. A constant entropy anchor
separates a complete commuting scalar gradient from the induction forcing.
`pair` returns the reduced electric field, potential, full field and entropy rate.
For the closed graph, dW_s/dt=J.T M E_s and the fixed-n logarithmic entropy
rate vanishes. For physical species, multiply the electron pair by beta and the
ion pair by −alpha. Their electric sum is (alpha grad pi−beta grad pe)/(ne).

This is an **explicit algebraic bridge/verification component**. Nodal thermals
remain lumped; they are not high-order physical thermal quadrature or actual
dual-volume means. The nonlinear face mean/product retains low-order pressure
error, despite the high-order field basis. The independently reconstructed
conservative density used by `VariableCompatibleCurrent` is not silently
replaced with this nodal density. Their full moving work/enthalpy coupling still
requires derivation and physical convergence checks. Neither class supplies
bulk advection, n_t thermal compression, a time integrator or vacuum boundaries.
