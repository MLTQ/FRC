# Parent-equation rate errors

Compares each conservative finite-volume rate against independent physical
cell averages of the parent species derivatives. Pressure rates differentiate
the declared numerical energy partition but compare with independent analytic
pe_t/pi_t; no equality is forced. Volume/dual, first-annulus, RMS and maximum
errors retain species discrepancies even when total E is conservative.

Manufactured forcing is explicit; wave derivatives are unforced linear modes
and omit quadratic terms by definition. Finite-amplitude controls are therefore
required. Wave error denominators use actual derivative perturbations, not a
large stationary background. Exact target arrays are archived for audit.
