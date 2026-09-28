# Independent curved, closed-flux equilibrium

`PeriodicEquilibrium` solves the continuum Grad–Shafranov equation with no
toroidal field and p(psi)=p0+lambda psi²/(2 mu0). Separation gives
psi=B0 R² F(r/R) cos(pi z/Z) and
F''-F'/s+(Lambda s²-k²)F=0. Its regular power series is computed without the MHD
finite-volume operator. A root solve sets F(1)=0 for the fundamental radial
mode. The fixed 48-term series is restricted to 0.5 <= Z/R <= 4; continuous ODE
residual and root accuracy are tested independently. This is a semi-analytic
reference, not a discrete equilibrium fitted to make the solver pass.

`fields` returns continuous SI psi, Br, Bz and pressure, including regular axis
limits. `state` samples exact corner flux, integrates continuum pressure over
annular cells (Gauss radial quadrature and analytic axial integral) and forms
total energy from the actual discrete B. Source and force cancellation are
then the responsibility of the MHD method; no equilibrium correction is added.

This fixture has curved closed poloidal flux, interior O-points and axial field
reversal, but repeats periodically inside a conducting cylinder filled with
fluid. It is deliberately initialized in equilibrium. It is not formation from
coils, a plasma–vacuum boundary, a stability proof, or confinement evidence.
Finite-time drift and initial force residual must be reported with refinement.
