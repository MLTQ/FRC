# Prescribed finite coils driving moving material plasma

Extends MaterialVacuumCurrent with analytic circular-loop fields evaluated at
the actual moving physical quadrature and a prescribed waveform I(t). The source
field remains separate from induced coefficients, so its exact curl is zero
inside plasma. A filament entering or touching the actual material polygon rejects. No
hardware circuit, exterior matter/conductivity or fixed plasma wall is added.

The total field in Hall and Lorentz forces is Bind+I*Bunit. Source material E is
-Idot*Aunit+I*v cross Bunit. Its actual current-metric load is subtracted from
the induced Faraday/Ohm forcing. The same full H,D and geometry derivatives
remain; physical Ohm residual restores the projected source-material E. Physical
laboratory E adds the exact -Idot*Aunit, keeping it distinct from its finite-space
projection. Current heat and pressure/entropy depend on actual J and species.

Source force is integral J cross (I*Bunit), including the r-weighted angular
torque. At fixed instantaneous state, the parent thermal/force calculation is
unchanged; only the mixed electromagnetic rows are solved again with the source.
This temporary repeated linear solve avoids changing an active archived parent's
implementation. It is not a split fluid/time integration. All stages use their
actual waveform time; unchanged1e-10 physical equation bounds remain enforced.

Reciprocal linkage lambda=integral J.Aunit dV=e*wcoeff.T*loadA. Its derivative
uses actual wdot, electric-map rate and the analytic moving-potential derivatives
d_r Atheta=Bz-Atheta/r, d_z Atheta=-Br. Therefore

    d(Efluid+Kcurrent+Wbind)/dt = -Idot*lambda,
    d(Efluid+Kcurrent+Wbind+I*lambda)/dt = I*lambdadot.

The prescribed ideal-filament self-energy baseline is excluded. Source work is
computed from waveform/linkage, never fitted to observed energy change. Both
budgets must be measured on complete coupled runs. Inherited upwind numerical
diffusion and measured current time sensitivity remain; no formation, physical
transport/edge or experimental validity follows from power cancellation alone.

Independent source groups now each contain fixed loop geometry and I_g(t).
`applied_source` supplies Bsrc=sum I_g*B_g, Asrc=sum I_g*A_g and the actual
partial time derivative sum Idot_g*A_g. All Hall/Lorentz/motional and electric
terms use those sums. lambda_g=integral J.A_g dV is retained for each group:

    d(Efluid+Kcurrent+Wbind)/dt = -sum Idot_g*lambda_g,
    d(Efluid+Kcurrent+Wbind+sum I_g*lambda_g)/dt = sum I_g*lambdadot_g.

All pure prescribed-source energies, including source/source mutual energy,
belong to the excluded source baseline. Plasma/source mutual work is explicit;
it is not obtained by fitting observed plasma energy. Single-group construction
and scalar observation fields remain compatible. For multiple groups, scalar
current/linkage fields are null and explicit group arrays must be used. Calling
`source` without selecting a group rejects ambiguity; `applied_source(time)`
returns the actual total source. The regular axis potential has its exact zero
limit. Exact matter-free vacuum and all original physical force terms remain.

The driven electromagnetic rows reuse the parent's `collision_electric`, which
is either the original eta*J or the paired Maxwellian force R_e/(e*n). Species
heat is already applied by the parent and is not added again. This keeps the
source-driven drag and heat from one instantaneous collision calculation; the
prescribed source is curl-free inside plasma and does not change J at fixed state.
The current isotropic collision closure has no explicit magnetic-field dependence.
