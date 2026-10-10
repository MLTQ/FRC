# Exact rejected radial256 arithmetic diagnosis

Rebuilds and hash-checks the original rejected minus input. Assembles unchanged
full physical forcing and analytic kernel separation, then performs the first
inner attempt at the original tolerance with an observation-only solver copy.
Saves complete residual vectors, proposed coordinates and the algebraic-path
difference prediction. It never returns or accepts a full plasma derivative.
The original complete five-attempt failure remains separately archived. The
observed first-attempt failure message/residuals must match that original
attempt exactly, and a separate copy tests input immutability.

This is not a tolerance relaxation or solver repair. Its purpose is to test
whether inconsistent floating-point association explains why balancing stalls.
All physical matrices and mass inverses are original; any remaining mismatch
is measured. Full fresh replay is required before using the diagnosis.
