# Reproducible wall-mode diagnosis

Archives the old rapid mode's independent physical reconstruction at quadrature
16/20 and Fourier/Bessel projection with axial/radial cutoffs. Reports captured
norms and continuum eigen-equation residuals; incomplete transverse capture and
coarse background errors prevent claiming a complete matched-mode proof.

`restricted` tests the saved original Jacobian in the exact magnetic/current K
metric, with/without the normal-fluid restriction. Controls are explicitly
incomplete boundary experiments, not time-evolution results. `calculate` forms
the actual complete StrongWallPlasma Jacobian independently through centered
RHS differences. Uniform rest and two contact difference amplitudes are used;
the full candidate is compared with the separately restricted original matrix.
The uniform maximum-real part is compared with its full oscillation frequency;
small finite-difference growth there is not interpreted as a physical unstable
mode. The original unconstrained probe's null-mode limitation remains separate.

Saves all physical fields, projection arrays, matrices, eigenvalues, leading
vectors, coordinate lifts, strict induction/trace maxima and complete immutable
source closure. Refuses overwrite. All arrays/non-time metrics require a
separate exact replay before claiming replay certification. No lab conclusion.
