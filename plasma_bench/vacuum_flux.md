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

`flux_components(..., magnetic_energy_flux="face")` selects the magnetic part
of the total-energy flux.
- **`"face"`:** the historical form. It uses the reconstructed face states:
  central (B²/μ0)v_n − (v·B)B_n/μ0, plus the LLF jump of B_t²/2μ0.
- **`"dissipative"`:** that LLF jump becomes the Poynting flux of the UCT
  dissipative field (`uct_electric(..., parts=True)`), paired like
  resistivity. In 1-D with a common speed the two forms are equal,
  ½a ΔB_t·B̄_t/μ0 = ½a Δ(B_t²/2μ0); in 2-D they differ by the
  corner-versus-face stencils.
- **`"ct"`:** all magnetic energy transport becomes E_θ×B_t/μ0. E_θ is the full
  CT field (ideal, dissipative, resistive and interface), averaged from the two
  face nodes. B_t is the mean of the two reconstructed tangential traces, or the
  wall/end trace on boundary faces.
  - The gas energy then changes by −∇·S_CT − dM/dt ≈ E·J, which is small where
    J ≈ 0. Near-vacuum gap cells no longer exchange large spurious energy with
    their gas.
  - Uniform states keep E = −v×B, so the flux and the jump conditions are
    unchanged.

The result carries `dissipative_electric` for the heating estimate.

**Fixtures** (`examples/dual_energy_verification.py`, `ct_flux` against
`entropy_reset`):

| Fixture | ct_flux | entropy_reset |
|---|---|---|
| Standing shock, ρ L1 | 0.39917 | 0.39917 |
| Standing shock, p L1 | 2.5859 | 2.5867 |
| Coil pulse, max difference from off | 0.0015 | — |
| Low-β adiabatic-invariant error | 0.0305 | 0.0305 |
| Fluid residual | 3e-16 J | — |
