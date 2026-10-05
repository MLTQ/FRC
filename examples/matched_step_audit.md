# Independent saved-stage substitution

Verifies source/raw identities and freshly evaluates every actual saved pilot
stage through the original complete physical equations. Reconstructs real Gauss
stage and endpoint equations independently from the saved coefficients and
requires the original stored residual bound. Retains full residual fingerprints
and induction checks. No nonlinear iteration algorithm or inverse proposal is
used to infer that a saved step passed. Accepted partial pilots can be audited,
with explicit distinction from completion of every requested step.

Warms the fixed inner-inverse cache with the same initial physical state as the
original pilot. All actual density-dependent operators stay physical. The audit
can be repeated to compare every numeric metric/fingerprint bitwise. It checks
equation satisfaction, not the pilot's continuum/time/mesh physical accuracy.
