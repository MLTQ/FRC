# Bounded upstream pressure trace

Purpose: solve the convex positive-k part of the actual current midpoint step
in `material_split_current`, retaining the curl metric's null gradients. This
is a numerical solver component, not a verification framework.

`box_minimum(A,b,limits,guess)` minimizes y.A.y/2-b.y with the stated symmetric
positive-semidefinite matrix and individual bounds. A bounded quasi-Newton solve
is followed by a working-set solve of the free equations. Newly hit bounds stay
fixed until the free equations settle; only then is one wrong-sign constraint
released. Null gradients remain null through a scaled least-squares free solve.
The function returns the pressure trace; it never clips plasma entropy or field.
Nonconvergence rejects instead of changing the force or adding a matrix penalty.

The earlier pulse stopped at step808 because a bound was repeatedly released
from its gradient sign while the coupled free Newton direction pointed outward
there. The allowed step length was exactly zero. Retaining a working set removes
that stagnation. Evidence is the actual failed-step matrix in
`outputs/material_split_failure_diagnosis`; the corrected1000/2000-step weak pulses now complete and show reduced temporal
current sensitivity. This does not establish strong-pulse or nonconvex-branch
validity; see docs/MATERIAL_EXPERIMENT.md.
