# Charge transport paired with the higher-order field metric

`CompatibleCharge` retains the closed periodic nodal graph and physical dual
volumes of `DualCurrent`, with natural Btheta(R)=0. Its stream is
S=(H_g g)/(dx dz), replacing the old diagonal-mass expression. Radial current
flux is S_previous_z−S and axial flux is the radial difference of S, padded by
zero axis/wall streams. Their graph divergence cancels locally.

These are **weak graph fluxes**, not claimed exact integrals of reconstructed
physical J through geometric dual faces. They satisfy
I_r=(M J)_r/dx and I_z=(M J)_z/dz and therefore give the same work against a
compatible scalar gradient as Ampere current and the full physical Hessian.
Using the old pointwise g stream with the new H would break this pairing.

`work` uses the full H_g metric. The old diagonal `magnetic_mass` attribute is
deliberately absent. Only zero wall Btheta is supported. A graph
pressure transport may use these fluxes to retain algebraic thermal work; this
alone does not raise its physical pressure accuracy or supply moving-density,
bulk enthalpy, stress or Poynting transport.
