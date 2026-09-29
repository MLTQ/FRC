# Physical finite-inertia current on higher-order compatible fields

`CompatibleCurrent` supplies constant-density physical electric mass M,
magnetic Hessian H and inertia I=mu M/(n e²). Ampere current solves
M J=−C.T H B and therefore pairs exactly with magnetic work. Physical electron
and supplied ion masses define mu and the Hall coefficient; no fitted mass or
current is used. The existing original-equation induction solver consumes this
contract without changing its residual tolerance or physical equations.

This component supplies instantaneous current/inertia only. It is not the
nonlinear moving finite-mass integrator, species transport or a fixed-density
substitute for formation. Constant n is an explicit linear verification domain.
Variable-density Hodge/inertia, pressure/work/thermal pairing and conservative
moving coupling must be derived and verified before production adoption.
The separate `variable_compatible_current` now verifies conservative variable-n
inertia and local energy derivatives. Its explicit graph pressure bridge still
has a nonconverging near-axis entropy-rate defect, so full moving adoption remains
open. This constant-density class and its archived linear studies are unchanged.
