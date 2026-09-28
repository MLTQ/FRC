# Independent magnetic shock verification

Verifies Rankine–Hugoniot mass, momentum, energy and induction fluxes, entropy
increase and the hydrodynamic zero-field limit independently of the numerical
operators. Differentiates the reference branches to check the cylindrical
momentum equation; checks Bernoulli/advected invariants and choking rejection.
The rejection radius lies below both branch-specific choking radii; the
lower-entropy upstream branch remains admissible farther inward than downstream.
Runs the actual bounded and vacuum-coupled 2-D solvers in the uncontaminated
interior window to require improving density/B errors and no false field sign
reversal with the declared numerical diffusion. Budgets remain independent
requirements. These small regressions supplement the larger recorded study;
they do not prove arbitrary-shock robustness or FRC formation accuracy.
