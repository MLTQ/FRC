# Primal momentum to regular dual velocity

Mass-weighted overlap of radial momentum/<r>, angular momentum/<r²> and axial
momentum, divided by overlapped rho, supplies nodal (vr/r, vtheta/r, vz).
The moments use exact annular <r> and <r²>. Constant regular compression,
rotation and translation are reproduced in the interior for constant density;
smooth variable profiles require independent refinement. No zero-radius
division or velocity floor is used. The boundary radial coefficient is set to
the explicit impermeable-wall value zero, not an inferred exterior flow.

This is a velocity representation, not another mass-flux calculation. The
actual continuity flux, including its declared numerical diffusion, separately
supplies dual n_t. Their finite-grid difference remains a numerical error to
measure. No density is independently evolved on the dual grid.
