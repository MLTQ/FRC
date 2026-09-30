# Experimental radial work pairing and its limits

The original LogFluid conservative pressure stencil has growing radial acoustic
modes even at uniform stationary density. Its finite-volume conservation alone
does not make its pressure and continuity derivatives an energy pair. The
preserved baseline source is 0aa417e588be92ab1ad5d1d5e1deeece7e6c2f02cb54ded16a01562849828fea.

## Physical work construction

At fixed n, physical kinetic energy has radial momentum Hessian

    H_r = M R_r^T W diag(1/n) R_r.

Here R_r reconstructs actual cell means of n v_r as regular r P(r²/2) fields,
and W includes cylindrical physical volume. B_r maps these same momentum means
to the existing normal particle flux r n v_r. Thermal transport supplies the
exact generalized radial face-work coefficient w_H for its actual n_t/s_t.

Let k_n and k_m be derivatives of the true global K with respect to actual cell
density and all three momentum means, including angular momentum. The radial
convective finite-volume flux is Gamma_r times the face specific momentum.
Transpose its divergence against k_m, keeping those specific momentum traces
fixed. Add the transposed continuity contribution against k_n. The resulting
w_K satisfies the actual radial convective kinetic-power identity. It includes
the changing-density kinetic term; omitting that term is not permissible.

The candidate source solves

    H_r a_r = −B_r^T (w_H + w_K).

It replaces the baseline radial pressure/pressure-hoop and centrifugal sources.
The conservative convective fluxes themselves remain. In the smooth continuum,
w_K corresponds to −M vtheta²/r times radial particle flux, so its negative
adjoint supplies the physical centrifugal source. A pressure-only preliminary
trial retained a different centrifugal discretization and failed rotating
balance badly; its source and measurements are preserved separately.

The actual radial-only kinetic-plus-thermal power then cancels to the iterative
solve/quadrature tolerance. Axial and angular momentum fluxes remain conservative.
The independent cell-energy account and all physical energy defects remain
visible. Axial work is not made exactly compatible, and the earlier uniform-
translation obstruction still applies. No post-step temperature/energy reset,
numerical diffusion, floor or artificial atmosphere is used.

## Linear and nonlinear verification

For constant n and entropy, define physical density/radial kinetic metrics
H_n/H_m and conservative divergence D. The baseline pressure map G gives
A=c[[0,−D],[G,0]]. Its diagnostic energy-adjoint counterpart is
G_ad=H_m^-1 D^T H_n. Direct spectra, central differences of the original full
nonlinear RHS, an independent matrix exponential and three RK4 step controls
separate spatial growth from time error. This linear control does not prove a
nonlinear force-accuracy or stability theorem.

The candidate physical radial metric is solved with a constant-density tensor
Gram preconditioner and checked relative residual. The larger radial acoustic
frequency also changes the explicit time-step restriction: the original nr8
contact step is too aggressive in the variable-density state. Half- and quarter-
step controls are retained rather than silently rerunning with an undisclosed
smaller step. Atomic stage rejection protects the last accepted physical state.

## Unresolved accuracy requirement

The independent solid-rotation residual improves sharply. However, a fully 2-D
rotating/compressing parent-rate comparison exposes a remaining radial defect.
Degree5 full acceleration relative L2 falls from .5015 at nr8 to .007346 at nr64,
but its FIRST-ANNULUS RADIAL component error rises toward about 6.7 times the
local physical radial acceleration. Axial and azimuthal component errors there
converge to about 8.54e-10 and 3.97e-12. The overall norm and shrinking axis-cell
volume must not be used to hide this local failure.

Independent force-term separation at nr64 gives first-annulus pressure and
kinetic-source relative errors6.44600/6.44732, versus2.15e-10 for conservative
convection. Thus both work sources need correction; the convective transport
already converges. Their compensating errors explain why good rotating balance
does not establish correct individual forces. See outputs/radial_force_terms.

The derivative of the reconstructed kinetic energy produces dual cell work
coefficients, whereas convective fluxes and the divergence-commuting particle
lift use different reconstruction maps. The observed error motivates a joint
analysis of those maps and their near-axis polynomial consistency; that is a
working diagnosis, not yet a proven root cause or correction. The candidate
must not replace the full moving plasma solver until this is resolved and
reverified. Field/current coupling, physical edge, transport and formation remain
open independently of the bounded contact/wave successes.

Candidate frozen source: 4f0c824a662e581dc78ddfbd81a2c9cc06342fb85d19d62fda44d0fb709ab3a9.
Raw rate/audit records: outputs/radial_work_rates. Spectra:
outputs/log_fluid_radial_diagnosis. Translation obstruction:
outputs/log_fluid_translation_probe. Original and revised local pilots and
time controls are separate directories; remote fourteen-case studies retain
their own platform provenance and must be audited on that platform.
