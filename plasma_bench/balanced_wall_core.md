# Balanced dual-residual wall refinement

Opt-in solver for the unchanged admissible positive induction equations.
The original operator uses preassembled (H T)y while physical reconstruction
uses H(Ty); mathematically equivalent actions can round differently. A correction
from only physical Faraday can reduce that defect while leaving the original
operator residual above its limit. When the latter fails, this solver uses the
arithmetic mean of the two dual residual correction vectors. Otherwise it uses
the original physical correction.

All original residual definitions and acceptance limits, inner CG tolerance,
three-refinement budget, physical inertia, admissible lift, and fixed inverse
remain unchanged. Neither residual is waived; failure is atomic. No equation,
state, rate or heat repair is applied. `refinement_trace` records every original
residual at each attempt for diagnostics. This is not a continuum accuracy fix.
The legacy core and all frozen source closures remain unchanged.
