# Common radial momentum: scope and equations

This is a separate smooth, unmagnetized fluid verification component. It has not
replaced the coupled moving plasma solver. The radial force defect in
RADIAL_WORK_NUMERICS motivated a change of state representation, not a fitted
force factor. Physical masses, pressure laws and reference inputs are unchanged.

## Why the old radial adjoint was inconsistent

The previous cell radial momentum used a regular reconstruction R, whereas its
continuity flux used a different lift L=QB. Even when both reproduced a smooth
polynomial exactly, solving `(Rᵀ W R) a = Lᵀ W f` did not reproduce that force.
The independent wall-compatible force `f=(r/Rwall)(1-r²/Rwall²)²` retains a
first-annulus relative error approaching 6.43 under degree5 refinement.

The map B from nr cell momenta to nr−1 interior radial flux faces has rank nr−1.
Replacing the cell kinetic metric with `Lᵀ W L` would therefore be singular.
The new state instead evolves the nr−1 independent interior faces themselves.
This is a new state space, not a relabeling or silent projection of old records.

## Physical state and work

The state consists of actual cell means n, mz=n vz, L=n r vtheta; interior
radial faces q=r n vr; two nodal physical log entropies; and an independent
conservative diagnostic energy account. Axis/wall q is exactly zero and z is
periodic. The radial physical particle momentum is Qq, where Q is precisely
the primitive lift used by continuity, divided by r. Thus radial velocity,
radial transport and radial kinetic work share one physical field.

With M=mi+me and actual reconstructed positive density n, the true kinetic
energy is `K=M |m|²/(2n)`. The radial Hessian is

`Hq(n)=Qᵀ W diag(1/n) Q`.

Axial/angular momenta retain conservative finite-volume transport. Their global
integrals telescope with the impermeable radial/periodic axial boundaries.
Continuity uses the actual q and centered axial momentum faces. The common
commuting lift obeys sampled `div Γ = −R n_t`. Both species use the existing
number-weighted log-entropy transport and physical pressure law; no temperature
is inferred or reset from the energy account.

The radial particle-momentum equation is

`M Hq q_t = M Qᵀ W/n [−Γr ∂r vr −Γz ∂z vr + vr R n_t + n vtheta²/r] − wH`.

Here wH is the sum of both species' actual thermal radial face-work coefficients.
The pressure contribution exactly cancels radial thermal work when paired with
the actual q. The centrifugal term is the original physical n vtheta²/r.
Physical quotient derivatives include cylindrical geometry; the radial momentum
derivative follows the same primitive's commuting derivative.

The physical radial Hessian is positive definite on the independent face DOFs.
The normalized iterative solve uses its physical Gram preconditioner and rejects
an unsuccessful relative residual. It does not fit density, current or mass.

## Conservation limitations remain explicit

A preserved intermediate prototype forced radial acceleration to cancel the
cell axial/angular convective kinetic-work defect. This gave only low-order
original-equation accuracy despite an exact radial total-work identity.
Its source, failed convergence test and measurements are archived under
`outputs/face_radial_paired_initial`. The diagnostic expression remains callable
but is not part of the evolved derivative.

The current physical-force method has exact radial pressure/thermal pairing,
but **not an exact nonlinear total-energy transport theorem**. Axial cell
momentum, its continuity lift, and radial strong advection still have different
discrete product rules. Physical total energy, the separately conservative
account, the RK4-integrated spatial discrepancy, and the remaining finite-step
chain-rule defect are all recorded independently. A smaller time step cannot
remove the spatial defect. No global energy correction or net-momentum force
is added to hide it.

Individual pressure and centrifugal accelerations in the manufactured reference
are nonzero at the outer wall; only their sum meets the fixed impermeable-wall
constraint. Their separate constrained projections therefore retain boundary
response errors against the unrestricted individual continuum fields. Raw
whole-domain norms are kept, while the original near-axis defect is assessed
separately. The independent polynomial force is individually wall compatible.

There is no added physical/artificial viscosity, resistivity, density floor,
energy repair or numerical atmosphere. Strong derivatives are cell-interior
derivatives of reconstructed fields; shocks/jumps have no conservation, entropy
or positivity theorem. Invalid reconstructed stages reject atomically. This
positive-density bulk component does not solve a heated free edge or replace
exact exterior vacuum with matter. Magnetic/current/Maxwell/drift and full
enthalpy/Poynting transport, physical closures and formation remain open.

## Reproducibility

`examples/face_radial_*` retain independent initial means/face traces, original
Euler rates, physical root-product/Legendre sampling, frozen source closure,
raw arrays, complete accepted histories, audits and figures. Long runs write
actual accepted-state checkpoints approximately every 5% of the requested
duration. Atomic metadata publication follows complete array serialization.
Resume validates source/raw hashes and bitwise initial state/rate on the
original platform. A test interrupts, serializes and resumes a nonlinear wave
and reproduces the uninterrupted arrays/history/energy ledger bitwise.
