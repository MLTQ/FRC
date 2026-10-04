# Species pressure and common charge-flux diagnostic

Uses frozen verified current-transport states and adds declared smooth,
nonuniform species entropy fixtures at the fluid points. It does not evolve
any state or implement a species solver. Independent physical reconstruction
measures pressure-current work, strong advective entropy work, their gradient
remainder, and total entropy drift at quadratures 16/20.

The continuum identity grad(p)/n=grad(gamma*h)-h*grad(s) gives the exact
pressure/thermal pairing only when current annihilates the enthalpy gradient.
The probe compares actual nonlinear gradients with compatible incidence
gradients of nodal analytic potentials. The latter use the existing curl-grad
identity; approximation errors are explicitly measured. It also checks the
analogous entropy-gradient annihilation. Neither substitution alone defines a
conservative coupled entropy flux or proves accurate parent-equation evolution.

The independent transport sampler is reused without modification. Each input
raw hash is checked; the entire original closure plus this probe is frozen in
a separate output directory. Existing records cannot be overwritten. No force,
heat or entropy repair is applied. Source fixtures and numerical integration
are fully specified in code; measured defects guide frc-7mm.4.
