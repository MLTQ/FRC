# Physical warm Hall field norms

Retains the existing independent fluid/magnetic wave norms, then compares
continuous Q1 electron pressure/density and reconstructed total-minus-electron
ion pressure/temperature with the full warm analytic fields. Norms use physical
annular volume and normalize by the perturbation, not background. SI RMS errors
are retained alongside relative errors. Nodal electron fundamental amplitude
and phase are fitted to independent cosine/quadrature-phase templates.

`pressure_difference` compares two Q1 thermal fields on nested meshes by
integrating over fine-cell subdivisions, normalized by the fine trajectory's
actual pressure change. Three-point quadrature exactly integrates these squared
bilinear polynomials. This is a continuous represented-pressure norm, distinct
from conservative fluid cell restriction or an exact continuum solution.
Callers must pass nested meshes with coincident boundaries and periodic z.
