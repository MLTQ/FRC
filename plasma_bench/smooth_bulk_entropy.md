# Smooth dual entropy transport

Uses the actual mapped particle flux Gamma to transport U by Gamma U_face/n_face.
Both conservative U and n are reconstructed from dual volume averages with
shared cubic interpolation. Reconstructing the ratio U/n as if its cell value
were its true average would introduce a nonlinear contact inconsistency.
The actual density derivative is still supplied solely by the primal mass flux.

Radial interpolation fits moments on the nonuniform dual x=r²/2 volumes;
axial interpolation is periodic. For uniform translation in z, the primal flux
map and interpolation commute, recovering the same high-order linear advection
of both U and dual n. Constant U/n follows any shared mass flux exactly.
These properties need independent continuum/contact tests, not just algebra.

This central smooth-reference flux has no upwind dissipation, limiter, entropy
inequality or positivity theorem. Nonpositive reconstructed n/U rejects without
floors, repairs or hidden fallback. Zero radial boundary flux and exact periodic
shapes are required. Full multidimensional transport still inherits primal/dual
representation and flux-map errors; fourth-order full physics is not claimed.
