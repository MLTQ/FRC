# Opt-in refined moving assembly

`RefinedRegularPlasma` reuses the complete frozen RegularPlasma RHS unchanged.
Only its operator factory chooses RefinedKernelInduction, with the same geometry
cache and fresh actual density at each stage. The original class, source closure
and studies are untouched. This separate entry point makes adoption explicit.

All original full residual checks, physical masses, common continuity, species
temperatures, boundaries and atomic failure contracts remain in force. First
attempt success should produce identical numerical arrays. A refinement fixes
linear-solver accuracy, not time-step stability, physical closures or vacuum.

The frozen refinement study verifies all 360 coordinate perturbations at two
amplitudes, compares to the separately derived physical linearization, and
replays both matrices bitwise. One complete coarse fast-wave trajectory is
identical to the original. Long nonlinear contact and general moving-state
conditioning are not yet verified; see REGULAR_PLASMA_VERIFICATION.
