# Conservative physical pressure, momentum and enthalpy flux

One shared face integrates M n v v+p I and axial angular momentum n r vtheta,
with M=mi+me and p=pe+pi from the physical log-entropy state. Normal moment
traces use the centered integer-index face reconstruction; transverse traces
preserve actual cell means. Log entropy is continuous on its original nodal
basis. Face pressures are evaluated from these density/entropy traces, not fit
from analytic pressure. Nonpositive density or pressure rejects without repair.

The finite-volume derivative telescopes in axial and angular momentum. Radial
and axial differences are formed separately before summation, avoiding loss of
a small radial difference when a large uniform axial flux is present.
Radial momentum includes (p+M n vtheta²)/r under physical quadrature. Axis/wall normal
particle velocity is exactly zero; radial pressure traction remains. No matter
is supplied outside the closed cylinder. Axial faces are periodic.

The energy flux is (K+sum gamma_s ps/(gamma_s-1))v. It advances an independent
conservative cell-energy account. The true reconstructed kinetic and two thermal
energies are evaluated separately; their local/global difference is not reset.
This flux is consistent and smooth but has no exact physical-energy/entropy
theorem with the current mixed reconstruction. No numerical dissipation is
added and there is no shock positivity theorem. Magnetic/current stress and
Poynting transport are not present in this hydrodynamic component.
The measured radial acoustic spectrum contains growing modes; this baseline
stencil is retained for diagnosis, not adopted as a stable full plasma solver.
