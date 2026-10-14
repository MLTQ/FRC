# Material dual-volume species entropy

Replaces the non-monotone skew entropy representation when requested by
MaterialThermalCurrent. The preserved high-order heated run developed electron
entropy below its initially uniform minimum, then failed the Coulomb-log branch.
This method is first-order upwind on explicit dual volumes. It introduces
numerical transport diffusion; no entropy or heat is clipped. Complete new
trajectories are required before an accuracy claim.

Each Q1 cell is split at its reference midlines. Quadrature particles belong to
their nearest corner; adjacent pieces form a material dual volume per node.
The selector S has one unit entry per particle. Actual fixed particles give
positive node masses m=S.T*dN. Entropy is s_quad=S*s_node, constant on each dual
volume. State representability is checked; no state is projected or reset.

The compatible gradient D has oriented entries (-1,+1)/h per electric edge.
Let q_edge=(Mn*w)_edge/h_edge. Weak Ampere implies I.T*q=0, with unscaled
incidence I. Species flux is c_s*q, c_e=-mi/M, c_i=me/M. Upstream transport is

    m*sdot = I.T [c_s*q*s_upstream].

This conservative spatial operator has a minimum principle and explicit graph
diffusion. Temporal monotonicity still depends on step size; maximum outflow/m
is returned. Pointwise normal current traces are not proved by graph balance.

For the actual EOS, hbar=S.T*(dN*h)/m is exactly dU/ds divided by node particles.
The matched pressure electric coefficient is c_s*s_upstream*D*hbar/e, its load
Mn times that coefficient. Thus e*w.T*F=sum m*hbar*sdot exactly. A constant
upstream anchor is split into the curl-free potential. Restoring
grad[c_s*(gamma*hbar-s*hbar)/e] recovers the continuum pressure field in the
smooth limit. Pressure work changes together with entropy flux; no energy repair.

Heat is integrated over each dual volume and divided by its actual dU/ds;
its physical work remains the original quadrature heat integral. No edge heat
is dropped. Both species share gamma and density, so temperature ordering within
a dual volume agrees with entropy ordering. Exchange has one sign throughout
that volume, conserves energy and produces nonnegative entropy. With positive
Joule shares, the initially common entropy minimum is a spatial lower bound.

The thermal representation and pressure forcing are lower order than degree-5
fields. Their spatial error, edge traces and numerical diffusion must be measured;
accurate integrated energy alone is insufficient.
