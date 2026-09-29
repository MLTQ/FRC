# Conservative current stress and full energy face powers

Adds mu/e times the charge flux transporting each nodal drift velocity
component to radial/axial momentum; angular momentum carries r wtheta.
The dual divergence-free charge flux is reconstructed normal-continuously and
integrated onto primal faces by scalar_current_flux. It is not the sampled
weak current asserted to be divergence free. The cylindrical radial source
includes the actual quadrature integral mu n wtheta²/r.

Total energy receives computed full E cross B and charge transport of
[mu v dot w - delta mu |w|²/2 + (alpha Hi-beta He)/n]/e.
Bulk advection of current kinetic energy is already in finite_mass_flux and is
not added twice. Scalar potentials remain in E even when their curl vanishes.
Each interior face is shared, with periodic endpoints and explicit zero normal
wall flow/current. Global conservative energy is a face-flux identity; local
stress/current-energy consistency still requires independent convergence tests.

Input radial/axial integrated flux arrays are local RHS work arrays and are
modified in place. The physical state and evaluated fields remain unchanged.
Nodal drift traces follow the existing current-thermal convention; finite-grid
interpolation discrepancy is explicit, not a physical collision/heating model.
