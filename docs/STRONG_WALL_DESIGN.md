# Candidate boundary-conforming bounded plasma space

The rapid stationary-contact mode has large physical outer-wall traces despite
the older weak natural Btheta(R)=0 and zero numerical fluid flux. This motivates
an explicit approximation space satisfying the *same intended continuum wall
conditions*. It is not an exterior/vacuum model or a fitted magnetic state.
Old solvers, data and all original checks remain frozen for comparison.

Let x=r²/2. At the outer wall Btheta=R w g Ez, with w the exact radial
histopolant trace. Define Cg=(w tensor I_z)g. The admissible magnetic space has
Cg=0, as well as the existing fixed axis/wall psi. Bulk q=n vr/r at the outer
wall is an inactive zero degree of freedom. Tangential velocities remain free.

Let H be the physical magnetic Hessian, F the original Faraday map and M the
electric mass. The H-orthogonal projection P onto the admissible B space gives
the variational curl Fw=P F. It preserves curl(gradient)=0 because F does, and
preserves div B exactly in the original poloidal representation; the new
constraint only restricts the toroidal component. This is the declared new
Galerkin spatial operator, not a tolerance relaxation or post-step state repair.

For the toroidal block, H=Hr tensor Hz /mu0. Thus the projection is
Pr=I-Hr^-1 w^T (w Hr^-1 w^T)^-1 w in radius, tensor the axial identity.
All matrices and trace constraints are explicitly saved/verified. For an
admissible B, P B=B and P^T H B=H B, so weak Ampere current
Aw B=-M^-1 Fw^T H B equals the old Ampere current. Its radial component at
the wall is zero because it has the same radial histopolant as Btheta and
the axial weak derivative commutes with this trace constraint.

The full inertial induction equations in the new space are
M E=I(n) Jt+f, M Jt=-Fw^T H Bt, Bt=Fw E.
Every declared equation retains the original strict residual tolerance. The
unprojected F E-Bt is separately reported as a boundary discretization defect;
it must not be called a passing old-incidence equation. Its H-dual annihilates
admissible magnetic tests. Equivalently, with A=-M^-1 F^T H and
K=H+A^T I(n) A, magnetic evolution satisfies the original K equation tested
only in the boundary-conforming space. This supplies the physical reaction
multiplier; it does not select a scalar correction to fix total energy.

Instantaneous physical magnetic/current energy is unchanged by the restriction:
B^T K Bt uses the same actual n and physical species masses. For an admissible
B, the boundary reaction has zero work. The constrained bulk wall q has zero
energy gradient wq=2 M x q/n, so its normal reaction also has zero work. Density,
entropy and other momentum rates are not repaired. Tests must verify actual
total power, not merely these formal identities.

Time integration must parameterize admissible magnetic and bulk states directly,
including all intermediate nonlinear stages. Projecting an already accepted
state to hide an error is prohibited. Analytic initialization has its own
declared boundary-conforming approximation and independently measured error.
Kernel separation must apply the new projected curl explicitly; the frozen
old difference-curl check must not silently be reused or disabled.

Required before adoption: original continuum rate/axis/wall/max convergence;
physical field and total-energy audits; matched stationary spectra and weak
reaction work; strict atomic failures; complete time/mesh-refined contacts and
waves. Unresolved interior current divergence/jumps and local enthalpy/Poynting
interfaces remain separate gates. The eventual moving vacuum edge and actual
coil/end conditions require their own physical boundary treatment.
