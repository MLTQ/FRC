# Numerical sixth-order magnetic dissipation

An opt-in candidate for the unresolved short-wavelength moving-contact growth.
This changes the spatial equations and is explicitly numerical, not a fitted
physical resistivity or transport closure. All original moving, Hall, inertial
and species terms remain. No stored state is filtered; all stages solve the
modified RHS at the original tolerances. It is not yet an accepted remedy.

Let F be the wall-admissible Faraday incidence, H the magnetic energy Hessian,
M the physical electric mass and A=-M^-1 F.T H the actual Ampere map. Define
K=-mu0 A F. K is nonnegative and M-self-adjoint. For actual J=A B, add

    numerical Ohm forcing = eta6 M K^2 J
    numerical heat Q6 = eta6 |P K J|^2.

Here P reconstructs electric vectors. eta6 has units ohm meter^5. This gives
J.T forcing = eta6 (KJ).T M KJ = integral Q6. Deposit heat into the existing
species entropy densities using S.T W Q6/(V h_s) and an explicit electron/ion
fraction. Native total energy is paired instantaneously; temporal drift and
continuously reconstructed energy are separate measurements. The high-order
adjoint heat deposition is not nodally positive; no clipping is allowed.

In the smooth, charge-free interior this is E6=eta6 Laplacian^2 J, giving
sixth-order magnetic diffusion. The nonnegative heat differs from pointwise
J dot E6 by a divergence flux; do not infer a physical local heating law. The
repeated discrete operator supplies variational higher-order wall conditions,
not measured material wall physics. Existing fixed psi, zero Btheta and zero
normal velocity, regular axis and periodic z remain. Exact vacuum and physical
coil/end boundaries are absent.

Coefficient choice belongs to the caller and must be recorded. Resolution
refinement must reduce eta6 to zero before interpreting collisionless physics.
`numerical_hyper_power` is distinct from physical Ohmic power. This candidate
does not establish that damping resolves the pressure/current inconsistency.
