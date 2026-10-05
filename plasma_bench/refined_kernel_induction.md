# Strict full-residual-driven inner solve

## Purpose
Candidate precision correction for the captured localized magnetic rejection.
The original solver controls an energy-weighted physical residual and then
checks a distinct full coefficient backward residual. The first can pass while
the second fails; neither contract is dropped or rescaled here.

## Contract
`RefinedKernelInduction.solve` preserves the identical compatible equations,
analytic kernel separation, original outer tolerance and all original residual
definitions. It first runs the unchanged InertialInduction solve. If any full
check fails, it repeats with a ten-times tighter inner tolerance, up to four
additional attempts. Only a rate passing every original outer check is returned.
No input is changed; exhaustion rejects atomically. There is no timestep retry,
new diffusion, fitted field, energy correction, altered mass or thresholded curl.

Geometry factors are shared by shallow copies; each attempt still evaluates the
current model's physical density-dependent inertia. Successful first attempts
must reproduce the original solver's arrays bitwise. Metadata separately records
inner tolerance and outer refinement count; descriptive conditioning ratios are
not substituted for checked residuals. Returned iteration count is for the final
attempt, not a total work estimate.

This is a new candidate module; frozen RegularPlasma still uses KernelInduction.
Captured-input, independent full coordinate linearization and one unchanged full
wave trajectory now pass; RefinedRegularPlasma is the explicit opt-in entry point.
Broader moving adoption remains open. A tighter solve cannot fix arbitrary floating-point
kernel recombination or provide long-time nonlinear stability by itself.
