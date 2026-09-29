# Matched physical thermal rate cases

Uses the previously verified conservative density, exact surface-flux continuity
and independently initialized compatible magnetic fields. Specific entropy is
sampled directly from the continuum reference at nodes; it is not fitted to
discrete current or adjusted to cancel a numerical force. Both species share
the actual physical density metric with their separate gammas.

`target` supplies smooth nonbarotropic, constant-specific-entropy barotropic,
and isothermal continuum pressure/rate targets. The latter two use zero supplied
bulk density rate so spurious pressure-driven magnetic generation is visible.
Uniform temperature is not assumed exactly representable by finite log-entropy
polynomials at variable reconstructed density; its initial temperature error
and magnetic battery must be measured. Constant specific entropy is represented
exactly. Reversed-current control changes the actual initial field and reference.

Protocol: degree3/5 × nr8/16/32/64, quadrature8/12 controls, n0=1e16/1e20,
reversed current, one barotropic and four isothermal cases: **18 rate cases**.
