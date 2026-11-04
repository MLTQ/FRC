# iec_field_check.py

## Purpose
Independent end-to-end analytic comparisons of the new tetrahedral electrostatic
solver, without unit tests or implementation-shaped assertions.

## Experiments
`python -m examples.iec_field_check NEW_OUT` runs 40/28/20/14/10 mm sphere and
40/28/20 mm spherical capacitor meshes. Uniform real positive volume charge in a grounded sphere has
phi=rho(R²-r²)/(6 epsilon0), E=rho r/(3 epsilon0) and field energy
2 pi rho² R^5/(45 epsilon0). A spherical capacitor has
phi=(1/r-1/b)/(1/a-1/b), E=rhat/[r²(1/a-1/b)] and C=4 pi epsilon0 ab/(b-a).
Nodal potential, volume-weighted field and absolute field-energy errors are
reported against those independent formulas, including geometric mesh error.

An electron released in the uniform-charge field has angular frequency
sqrt(e rho/(3 epsilon0 m_e)). Full 3D drifts through the actual mesh over one
analytic period are compared with the analytic endpoint and analytic mechanical
energy, at two timesteps per spatial mesh. No particle position/energy repair.
This fixture holds the analytic charge fixed; it is not a physical plasma or
a claim of self-consistent PIC convergence. PIC budgets and particle refinement
must be measured on actual dynamic experiments separately.
