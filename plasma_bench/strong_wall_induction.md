# Strict induction in the wall-conforming space

Solves M E=I Jt+f, M Jt=-Fw^T H Bt, Bt=Fw E with StrongWallGeometry's explicit
Fw=P F. Uses the same positive magnetic-space iteration and fixed tolerance,
with inner refinement until all declared physical/coefficient/backward residuals
pass. Exact curl-kernel fields are split analytically and restored in full E.
The difference-first curl is explicitly projected by P; the old frozen raw
incidence solver is not silently changed or called a passing check here.

`unprojected_curl_relative` is a separately reported discretization defect,
excluded from the new Galerkin equation residual and never asserted zero.
Independent original-equation convergence and boundary-work checks are mandatory.
No altered mass, state correction, tolerance relaxation or energy/heat repair.

The initial larger-space Krylov pilot rejected a small inadmissible toroidal
rate trace during a contact mixed-direction test. Its archive is retained.
WallInductionCore now parameterizes every Krylov rate in the exact trace null
space; complete unreduced Fw/Ohm/Ampere bounds still apply. This is an internal
linear coordinate change, not projection of an accepted derivative. A completed
pre-change contact and its bitwise replay remain separate frozen evidence.
