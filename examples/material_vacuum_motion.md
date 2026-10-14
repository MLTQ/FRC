# Magnetic exterior under actual nonaffine boundary motion

Runs the new shape-dependent magnetic exterior through a complete prescribed
material-boundary trajectory. This is the exterior component needed for a free
plasma, not a plasma evolution or FRC formation run. Vacuum has no matter or
conductivity. The inner boundary expands25% radially at the midplane and15%
axially at the poles; corners stay fixed. All exterior geometry is recomputed.

A fixed physical point dipole inside the excluded plasma domain supplies exact
boundary psi and independent vacuum Br/Bz. A separate64-point integration of
the explicit analytic field on the actual piecewise-linear plasma boundary gives
exterior energy by Green's identity, Maxwell traction and the moving-domain
energy rate. None of these references is used to fit the vacuum solution/force.

The derivative at held material flux and the rate with the actual changing
dipole boundary flux are distinct. Both are recorded: geometry work comes from
the analytic energy shape derivative; total rate adds the boundary flux work.
The independent fixed-field rate is minus magnetic-energy density integrated
against the moving boundary's normal velocity. This exercises the exterior
force and the Hessian derivative required by subsequent finite-current coupling.

All states, fields, forces, matrices and frozen source are saved. Existing
provenance (`source_sha256`)/resource guards are reused; no unit tests or new verification
framework. Field, traction, energy and power errors measure physical accuracy;
agreement of two algebraic power evaluations alone is insufficient. Parent
fixed-cylinder multipoles and exact polynomial outer matching remain in use.
