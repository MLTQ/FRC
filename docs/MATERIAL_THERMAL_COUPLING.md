# Moving current, separate species entropy and heat (2026-10-14)

This extends the material current solver from d6983b8f with actual electron/ion
entropy transport, pressure work in Ohm's law, prescribed scalar resistance,
Joule heat and conditional Coulomb temperature exchange. Two complete heated
expansions finish. It remains the axisymmetric toroidal-field/meridional-current
sector, with exact empty exterior and freely moving sides/ends. It is not yet
poloidal FRC formation, moving coil/vacuum coupling or a validated heated edge.

## Equations and assumptions

Use M=mi+me, mu=mi*me/M, alpha=me/M, beta=mi/M, w=vi-ve. Species order is
electron, ion, and their relative velocities are c_e*w=-beta*w and c_i*w=alpha*w.
With s=ln(T[eV])-(gamma-1)ln(n), h=e*T/(gamma-1), gamma=5/3,

    D_v s_s = -c_s w.grad(s_s) + Q_s/(n*h_s),
    E_pressure = [alpha grad(p_i)-beta grad(p_e)]/(e*n),
    E_resistive = eta*J,  Q_J=eta*|J|².

All previous finite-mass Hall, changing-density/current inertia, vorticity,
magnetic force and current-stress terms remain. Material pressure uses the
actual changing density and both physical temperatures at every RK4 stage.
Constant-entropy restrictions are removed only in the new thermal subclass,
where the corresponding transport/pressure terms are supplied explicitly.

The default friction-heat split is beta to electrons and alpha to ions, a
declared mass-weighted partition; an explicit alternative fraction is supported.
The selected eta=0.001 ohm-m is constant input, not an inferred Coulomb tensor.
Optional temperature exchange reuses the existing finite-mass Maxwellian Landau
coefficient: Qe=1.5*n*e*nu_E*(Ti-Te), Qi=-Qe. Its assumptions and primary source
are in [coulomb_exchange_wall_plasma.md](../plasma_bench/coulomb_exchange_wall_plasma.md).
Unsupported branches/nonpositive Coulomb logs reject without clipping.

The first implementation adapts the former skew pressure/current pair to material particle weights:
P=(N.T*dN*N)^-1*N.T*dN, A*s=dN*w.grad(Ps), K=(A-A.T)/2,
sdot=-c_s*K*s/dN. Its electric pressure load is the exact work adjoint,

    F_s = -c_s/(2e) integral dN G_e.T [h*grad(Ps)-s*grad(Ph)],
    e*w_coeff.T*F_s = sum dN*h*sdot_s.

Thus thermal and field/current work cancel without a global repair. Compatible
gradient charge work gives spatial conservation of species entropy for the
source-free transport. Constant anchors remove exact gradients from induction;
the thermal and self-drift scalar potentials are restored in the physical
electric diagnostic. The latter also subtracts v x B to return laboratory E.
Entropy itself remains a quadrature state; it is not reset to a projected profile.

Local heating is Q/(n*h), so its internal-energy work equals Q at actual material
quadrature. No heat is suppressed near the edge. No floors, atmosphere, changed
masses or evolved-state energy corrections are introduced. The original skew
option has no added diffusion; its upwind replacement below introduces explicit
transport diffusion to address the measured entropy undershoot.
The same 1e-10 physical Faraday/Ohm residual bound remains. Linear refinement
solves the original complete mixed equations, including pressure and resistance.

## Complete physical trajectories and time comparison

Both use the compact n,T,B,J initial profile documented in
[MATERIAL_CURRENT_COUPLING.md](MATERIAL_CURRENT_COUPLING.md): R=35 mm,
half-length=70 mm, nc=1e18/m³, Ar+/physical electron masses, initial equal species
temperature profile, zero initial motion. The 4x8 fluid/field mesh uses the
existing degree-3 representation and 6-point quadrature. Motion and fields are
fully evolved for 7.121769 us. This initial current still has the previously
measured 32.2% spatial error; this time comparison does not eliminate it.

| Measured quantity | 1000 steps | 2000 steps |
|---|---:|---:|
| Maximum relative total-energy change | 2.3162e-10 | 1.4247e-11 |
| Maximum independent scalar-virial error | 3.3230e-9 | 8.1381e-10 |
| Observed Joule heat, nJ | 97.2740 | 97.2304 |
| Final particle-mean Te, eV | 0.370890291636 | 0.370890291508 |
| Final particle-mean Ti, eV | 0.369006374982 | 0.369006375011 |
| Wall time, seconds | 71.1 | 142.4 |

Halving dt reduces the total-energy error 16.26-fold. Final B/J relative L2
differences are 2.16e-7/3.97e-7. Electron/ion particle-weighted temperature
differences are 5.63e-10/1.39e-10 relative, maximum absolute 3.49e-9/2.07e-10 eV.
Final position differences are at most1.42e-11 m. These are time comparisons
at one coarse spatial representation, not a spatial accuracy certificate.

The independent continuum virial remains I''=4E-2Wb with both species internal
energies and current kinetic energy included. Its separate trapezoidal magnetic
integral explains its approximately fourfold time-error reduction. The Joule
and exchange integrals also use trapezoids; their convergence is distinct from
the actual RK4 species states. No integrated budget is fed back into evolution.

For 2000 steps, final Wb=2.15127e-12 J, Kj=3.80511e-11 J, peak B=0.158508 uT,
peak J=29.6013 A/m². The observed exchange transfers0.377111 nJ from electrons
to ions over the trajectory. Both species cool overall as the plasma expands,
while friction produces a warmer electron population. Maximum physical row
residuals are 1.87e-13 (Faraday),1.95e-16 (Ohm); maximum pair-power imbalance
8.17e-15 W and weak charge residual2.15e-16. Sampled total entropy production
is positive. None of these establishes a physical pointwise boundary closure.

## Accuracy and regime limits

The independent initial B/J projection comparison remains essential. Increasing
the existing representation to degree5 with10-point quadrature on8x16 reduces
initial B/J errors from5.21%/16.15% to0.0683%/0.2096%, with drift error0.3642%.
Magnetic/current energy errors are +1.052e-5/-2.520e-5 relative; independent
initial acceleration error is0.0006707 of C/R. This is initial-data evidence only.
The full nonresistive high-order trajectory is still running at this checkpoint;
do not claim its completion. The first high-order heated attempt failed below.

The completed heated runs encounter Coulomb logs as low as2.340 during evolution;
the final maximum-drift sample has log2.782 and drift/electron thermal speed8.12.
The helper's Maxwellian weak-coupling estimates are not reliable closure evidence
there. Prescribed eta plus conditional temperature exchange is not full physical
transport. General kinetic/neutral/ionization behavior, pointwise current
continuity, normal relative-drift traces and heated-edge regularity remain open.
Poloidal magnetic fields, moving vacuum and finite coils remain absent here.

## Preserved high-order failure and transport correction

The6x12,degree5,q10 heated skew run failed at attempted step187 (186 accepted,
1.32465us). Its frozen source and failure_stage_analysis.json are preserved.
The initial minimum Coulomb log was1.0261. At the last accepted state it was
0.01581; the next RK stage gave -0.04642, so the unchanged physical branch
rejected. At that corner, n remained8.20e10/m³ while Te fell from1.887e-5 eV
initially to9.623e-6 eV and drift reached about172km/s. Electron entropy fell
from the initially uniform -27.6310 to -28.3046, then -28.3461 at the rejected
stage. This is a numerical undershoot, not merely a physical Coulomb limitation:
the continuum species advection, nonnegative Joule shares and thermal exchange
preserve the initially common global entropy lower bound. At a minimum across
both species, relative transport cannot reduce it and exchange heats the cooler
species. The old skew operator lacks that minimum principle.

MaterialDualEntropy therefore introduces explicitly **first-order upwind
numerical transport diffusion** on positive material dual volumes. The quadrature
particles are assigned to the nearest reference corner, and entropy is constant
on each dual volume. Compatible Ampere covectors give an oriented, divergence-
free nodal particle-flux graph. With m=S.T*dN and q=(Mn*w)/h_edge,

    m*sdot_s = I.T [c_s*q*s_upstream].

For the actual EOS, hbar=S.T*(dN*h)/m. The matching pressure load is
Mn*(c_s*s_upstream*D*hbar/e), so its electric work exactly matches the thermal
transport work. The scalar pressure potential is restored separately. Actual
Joule/exchange heat is integrated over each dual volume and divided by dU/ds;
no edge heat, collision term or Coulomb-log restriction is removed. Since both
species share density and gamma, exchange has the same sign throughout each
dual volume and increases physical entropy. No state is clipped or repaired.

This spatial operator has the required lower-bound property, subject to temporal
resolution. The actual outflow/m and minimum species entropy are now recorded.
It is lower order than the degree-5 fields and needs spatial accuracy assessment;
energy conservation alone will not establish it. The complete replacement run
is in progress. The old skew runs remain explicitly selectable and archived.
The replacement has passed the failed run's1.32465us point: at1.79469us (step252),
minimum electron/ion entropies are -27.6304585/-27.63102098 versus initial
-27.63102112, minimum Coulomb log1.02710, and transport outflow*dt is0.05447.
These are partial observations; completion and time/spatial accuracy remain
unproven. At that point total energy change is1.73e-7, larger than in the coarse
skew comparison, so the new method needs its own temporal assessment.

## Evidence and live continuation

The [inspected updated figure](../outputs/material_thermal_comparison/material_current_expansion.png)
shows actual species temperatures, heat, field/current energy and final B/J maps.
Its PDF, `comparison.json` and `time_comparison.json` retain measured quantities.
The final B color scale uses the actual decayed field, not its initial amplitude.
The earlier comparison_r4 image remains as the inspected initial-scale version.

Completed thermal sessions96593 and84105 exit0 were collected. Plot observers
95847 and26949, and final-state comparison91296 also exit0. Completed raw sources
and states are in outputs/material_thermal_free_r4_t{1000,2000}_eta1e3.
Their frozen source closure is fea72a1041a501f76842a5c15a912ed66f8f1d3505e810435b697a9c78297e3d.

Failure26494 is terminal1 with exit collected; failure-stage observer99229
is terminal0. LIVE science handles:43078 (nonresistive8x16,p5,q10,2000 steps)
and60494 (heated6x12,p5,q10,1000 steps with upwind transport). Re-poll these exact
handles and inspect terminal summaries before a restart or source edit.
Do not refactor the base geometry or thermal source closure while either uses it. The new shared geometry
temporarily duplicates the base implementation to preserve its running source.
After completion, consolidate it and compare complete spatial trajectories;
then advance poloidal/moving-vacuum/finite-coil physics and physical transport.
The full FRC goal and frc-vacuum.3 remain open. No unit tests or new verification
framework were added or run. The older interactive3D viewer is unchanged.
