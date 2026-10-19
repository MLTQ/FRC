# Conservative gas-chamber MHD fluxes

`flux_components` uses shared limited primitive reconstruction, Rusanov fluid fluxes,
and four-state LLF upwind constrained transport. Mass, both momentum components
and total fluid-cell MHD energy advance by flux differences plus the cylindrical
radial total-pressure source. The induction flux has the same one-dimensional
Rusanov limit. No separately tuned compressive magnetic diffusion is added.

Two axial reflected trace layers prevent periodic end communication. Normal
velocity reflects at gas boundaries; magnetic field does not reverse at a
magnetically transparent radial wall. Boundary normal B comes from the common
staggered field. Tangential extrapolation/reflection supplies the slip-wall
trace. The ghosts are only reconstructed boundary data, not fluid/vacuum states
being evolved. Only the symmetry axis enforces zero corner electric field.
The actual boundary fluid flux uses exterior tangential B from `vacuum_trace`
for Maxwell traction and energy transport; simply reflecting interior B would
omit interface-current force/work. The induction boundary still uses its shared
normal B and tangential slip velocity.
LLF interface electric dissipation is included in induction and Poynting flux
on **all adjacent fluid faces**, including tangential transport along the
boundary layer. Restricting it to outward boundary flux misses internal energy
redistribution. It is not counted as physical eta J heating.
`spatial_order=1` supplies piecewise-constant conservative fluxes and CT electric
fields as an admissibility baseline. Second-order remains the preferred update;
any blending must apply to all fluxes, induction and metric sources together.
It returns `FluxUpdate` with shared face fluxes, corner electric fields and the
radial metric source. `flux_rhs` retains the assembled derivative API by calling
its `rhs`; local limiting selects shared entities before this same assembly.

Physical eta J at corners comes from the full plasma/vacuum stiffness and its
actual plasma-node areas. Its Poynting flux enters total energy, with no separate
Joule source or energy correction. Returned outward total-energy flux supplies
an independent fluid-domain conservation diagnostic. Global response energy
also depends on the exterior field and is checked separately, not forced exact.

`flux_components(..., local_bounds=False)` forwards the optional neighbour-minimum
trace bounds to `Reconstruction`.

`open_ends=True` sets the axial gas boundary per radial row:
- **Outflow rows** (end-cell velocity pointing out of the domain) use
  zero-gradient ghosts.
- **Inflow rows** use the exact mirrored reflecting ghosts, so nothing flows in.
- The magnetic boundary is unchanged.

Verification (`examples/open_ends_verification.py`): a uniform magnetized flow
at 0.5 stays exactly uniform downstream of the upstream rarefaction head, and
mass and energy close exactly with outflow counted. The total outflow is 1–1.5%
above the downstream-only value, because the stalled blocked end
back-flows slightly.
