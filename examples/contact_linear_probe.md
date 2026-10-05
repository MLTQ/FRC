# Stationary-contact spatial diagnostic

Removes the uniform axial motion from the exact translating contact. The
continuum state has zero velocity/current, constant total pressure and axial B,
and a nonzero curl-free electric field. This is a distinct stationary analytic
verification fixture, not a modified accepted trajectory. A uniform rest state
with the same mean density/pressures is the control.

On a deliberately small 360-coordinate mesh, centered differences of every
actual strict full RHS coordinate form the frozen spatial Jacobian. Records
base stationarity, all eigenvalues, a leading vector, eigensolve backward error,
mixed-direction linearization check and unchanged induction maxima. Two contact
perturbation amplitudes check finite-difference dependence. Frozen raw matrices
and source closure are saved; no time integrator enters this diagnostic.

Growing eigenvalues would show a spatial linearization issue (or continuum
instability represented by that discretization), independently of a chosen
time method. They alone do not distinguish those two possibilities, nor prove
the moving nonlinear trajectory's cause. An independent continuum linearized
operator, resolution/mode controls and perturbation evolution would still be
required. The literature on density-wave local linear instability is motivation
only: https://arxiv.org/abs/2009.13139 is not evidence about this plasma operator.
