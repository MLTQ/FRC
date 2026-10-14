# Material species-relative entropy and pressure work

Purpose: carry electron/ion entropy with their own velocities on the moving
material fluid, pairing that transport with the pressure term in Ohm's law.
This adapts the existing RegularSpecies physical skew pair to fixed material
particle weights and finite, freely moving ends. It is not a new verifier.
MaterialThermalCurrent integrates it with actual current/fluid/heat evolution.
Two complete heated runs and measured limits are in MATERIAL_THERMAL_COUPLING.md.

Species order is electron, ion. Their velocities relative to the bulk are
c_s*w with c_e=-mi/(mi+me), c_i=me/(mi+me), w=vi-ve. Physical specific entropy
is kB*s/(gamma-1), with s=ln(T[eV])-(gamma-1)ln(n). Internal energy per particle
is h=e*T/(gamma-1). Entropy remains an actual quadrature state, not a prescribed
temperature or a reset to a projected profile.

N is the existing compatible scalar interpolation and dN the fixed particle
weights. P=(N.T*dN*N)^-1*N.T*dN projects a quadrature scalar only for computing
its gradient. G is the physical compatible gradient on the evolving geometry.
For a given actual drift, A*s=dN*w.G(P*s) and K=(A-A.T)/2. The evolution is
sdot_s=-c_s*K*s/dN. For weak Ampere current, A*1=A.T*1=0 up to roundoff.
Thus constant entropy, its particle integral and the frozen-current quadratic
entropy integral are preserved spatially. This does not establish local
pointwise current continuity, interface transport accuracy or shock stability.

The pressure load is the exact work adjoint of this transport:
F_s=-c_s/(2e)*integral dN G_e.T*[h*grad(Ps)-s*grad(Ph)]. Therefore
e*w_coeff.T*F_s=sum dN*h*sdot_s, which cancels the corresponding field/current
power. Constant scalar anchors split exact gradients away from the reduced
induction forcing, following RegularSpecies. The full electric potential is
c_s/e*grad(P[gamma*h-s*h/2-ha*s/2+sa*h/2]); adding it reconstructs the pressure
field. In the continuum the sum is c_s*grad(p_s)/(e*n).

No entropy projection is imposed on the evolving state. No exterior matter,
heat suppression, numerical diffusion, floors or energy correction is present.
Separate physical collision/heating sources must be coupled explicitly before
claiming a thermal plasma model; this operator alone contains neither.

`pair` requires a current geometry, fluid state and drift from the same stage.
It returns the actual entropy rate, reduced pressure load, full scalar-potential
electric contribution and physical work/entropy rates. `project`, `gradient`
and `gradient_transpose` implement one matched physical projection and adjoint.

This skew option is preserved: the high-order heated6x12 run developed an entropy
undershoot and rejected at step187 on a nonpositive Coulomb log. MaterialDualEntropy
provides explicit upwind transport instead. `heating_rate` retains the old local
point heat and `potential_gradient` exposes its existing projection to composition.
