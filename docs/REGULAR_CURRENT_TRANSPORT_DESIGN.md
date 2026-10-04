# Historical current-transport derivation

The candidate below has now been implemented as a separate extension. See
CURRENT_TRANSPORT_NUMERICS.md and CURRENT_TRANSPORT_VERIFICATION.md for its
actual discrete contract and measured scope. The derivation below is retained
as the pre-implementation checkpoint.

The regular coupling checkpoint verifies actual fields and force/work adjoints.
The following continuous weak identity suggests a compatible next implementation;
it is not a verified discrete transport operator or an adoption claim.

For physical current J, a current test field t, physical density n and velocity
v, define Dv as the full cylindrical covariant directional derivative. With
axisymmetry, its transverse components include −vtheta*Jtheta/r and
+vtheta*Jr/r. Their axis limits must follow the regular fields, without a floor.
An advection bilinear form is

    <t,Fv> = mu/(2e²) integral [t·Dv J − J·Dv t]/n dV
           + mu/(2e²) integral (t·J)*n_t/n² dV
           + mu/(2e²) integral_boundary (v·normal)*(t·J)/n dS.

At t=J the skew volume part vanishes. Its remaining volume work is precisely
the density exchange exposed by the present component; its boundary work is
the current kinetic-energy advection flux. Together with the existing full
changing-density Ohm term, this is the desired temporal work cancellation.

Integration by parts gives the actual strong field

    Fv = mu/e² [Dv J/n − J*(Dv n)/n²]
         + mu*J*[n_t+div(nv)]/(2e² n²).

The first bracket is mu/e*(v·grad)w for w=J/(en). The second vanishes for
the original continuity equation. For discrete reconstructed fields it is a
local consistency defect whose convergence must be measured, not a fitted
global cancellation or proof of original-force accuracy. Actual continuity
must remain the common regular-fluid derivative.

The matching bulk drift force comes from −div(mu*n*w*w). Its weak work is
integral mu*n*ww:grad(v) minus boundary traction work; the Ohm stretching term
mu/e*(w·grad)v has exactly the opposite current-energy work. The actual velocity
map, its covariant derivatives and their volume adjoints are required, including
axial density dependence and dependent transverse axis traces. Symmetric stress
must preserve actual angular momentum and axial force as well as energy.

Implementation gates remain substantial: shared element-interface fluxes,
the impermeable wall's numerical versus reconstructed trace, normal current
and drift traction at boundaries, the self-drift term, and exact species
pressure/entropy work. A pointwise physical interpolant need not satisfy strong
continuity or div(J)=0 exactly just because compatible weak identities hold.
Do not drop these defects or boundary powers. Test independent original
advection/stretching/divergence rates separately near the axis and wall, then
the complete parent equations and unforced trajectories. Existing frozen source
closures must remain replayable; extensions should be separate until verified.
