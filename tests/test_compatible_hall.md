# Compatible Hall and coupled fluid verification

Independent physical differentiation verifies that coefficient Faraday curl
commutes with the reconstructed E/B fields, including the first annulus.
The original nonlinear snapshot now checks actual first-annulus convergence,
not only shrinking-volume norms. True physical quadrature, sparse Hessian and
the separately implemented energy gradient agree under the declared convention.

Nonlinear implicit midpoint is checked against independent adaptive DOP853
for temporal order, and against quadratic magnetic work, positive resistive
heat and conservative Poynting transport. Failed iterations preserve inputs.
Full HallMHD steps advance compressible flow and magnetic fields with positive
pressure and conserved mass/energy; a failed magnetic trial leaves every
accepted physical state/time/heat/iteration total untouched.

Explicit energy metadata is tested through construction and propagation.
Zero Hall/zero resistance matches the ideal RHS and two ideal half steps
exactly (not a claim of bitwise equality to one SSPRK2 full step). Continuum
wave initialization retains true conservative averages and quadrature controls.
Actual finite-pressure fast Hall propagation is refined on three meshes at a
separate common time cap. Invalid axis/wall/periodic data also fail before a
zero-transport shortcut can discard them. Full branch propagation, nonlinear
mesh/time studies and physical formation
remain separate measured requirements beyond these focused tests.
