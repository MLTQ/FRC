# Joint finite-mass current midpoint

Unknown magnetic/current state and both thermal entropies share one nonlinear
solve. Midpoint magnetic/current variables preserve quadratic B/current energy;
temporal logarithmic means preserve each pressure transport entropy identity.
The Hall matrix uses delta B*, including convective inertia. Pressure includes
both species with their exact mass fractions. Physical source equations are

    M dj_mid + dt/2 C.T H C E = 0
    M E - 2/dt M_inertia dj_mid = X(delta B*) (jold+dj_mid) + M Ep.

Solving for current increments avoids subtracting a large initial inertial
source. Both Picard transport and its residual use `dj=current(db/2)` with
`db=dt C E`, `jmid=jold+dj`, and inertia `2 M_inertia dj/dt`. Applying the
linear current map directly to the increment avoids catastrophic cancellation
in `current(bmid)-current(bold)`. No term or small current is dropped, and the
Ohm residual remains relative to its physical terms at the requested tolerance.
The independently recomputed increment current also checks the mixed solve.

Electric/current mass scaling and the existing checked sparse solve are reused;
their numerical cache may change on rejection. A separate stored-state check
uses the actual endpoint average and normalizes Ampere error by the termwise
absolute matrix-product envelope `|M||jmid|+|C.T||H||bmid|`. Faraday storage error
is normalized by `|bnew|+|bold|+|db|`. Their maximum is returned as
`stored_state_residual` and must pass the same tolerance. These backward errors
expose representational roundoff under a large bias; they do not replace the
small-signal Ohm test or imply that perturbations below storage precision resolve.
Both thermal equations are checked independently. Full electric and
conservative face energy include restored scalar potentials. No energy reset,
physical diffusivity adjustment or hidden retry is used. Inputs are copied;
invalid geometry, nonpositive states or residual failures reject the trial.
Returned radial/axial face powers are integrated watts, with explicit boundary
faces; they expose the actual flux used in the conservative cell increment for
independent continuum comparison and archival replay.

This is a constrained v=0/fixed-n component with exact physical electron mass
(or an explicit consistent massless control), not a complete moving finite-mass
solver or a proof of unconditional stability. Actual convergence, physical
current-energy/entropy replay, limits, local flux defects and failed fixtures
must be measured before integration into moving mass/momentum/energy.
