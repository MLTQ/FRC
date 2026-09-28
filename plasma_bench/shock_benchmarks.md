# Local standing shock in the actual axisymmetric kernels

The exact cylindrical steady branches are embedded in a 2-D unit chamber.
`regularized_profile` keeps that solution only on 0.38<=r<=0.85 m, clips smooth
branch data beyond it and makes radial velocity regular at the axis. Velocity
smoothly tapers to zero between r=0.85 m and the wall, avoiding an unrelated
vacuum-producing wall rarefaction from the original inward supersonic inflow.
`wall_taper=False` retains that remote stress configuration for a contamination
check; it is not used to claim free-vacuum robustness. These
outer parts are declared verification initial data, not part of the analytic
solution. Compare only a central annulus and axial window before disturbances
from those extensions or chamber boundaries reach it.

`initial_fields` integrates gas conserved variables and r Bz with independent
Gauss quadrature; the initial shock must align with a radial face. `standing_shock`
constructs either the conservative bounded reference or the actual split solver
with infinite magnetic exterior. The latter harmonically extends the specified
plasma flux into vacuum. End/interface surface-current effects are outside the
short-time interior comparison. Extending time/windows requires fresh boundary
contamination checks, not an assumption that the entire chamber stays steady.
The coupled factory exposes `midpoint_method` so the unchanged implicit physics
can be compared between direct and cached-vacuum fixed-point iterations.
The explicit `shock_diffusion` option compares the original centered reference
with vanishing numerical resistivity. Its heat is not counted as physical Ohmic
heating; the shock's independent jump/steady branches remain unchanged.
`inner` changes the nonanalytic inner extension, and `half_length=2` doubles the
distance to reflecting gas ends at unchanged dr=dz. A correspondingly larger
matching rectangle still encloses the source sphere. These are explicit boundary
contamination controls, not alternative exact global steady solutions.
