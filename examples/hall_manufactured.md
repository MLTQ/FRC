# Nonlinear Hall manufactured snapshot

Prescribes smooth finite-amplitude axisymmetric fields and positive variable
carrier density on a bounded cylinder. psi=B0 r²/2+A r² F cos(kz),
Btheta=C r F cos(ell kz+.37), n=n0[1+delta F cos(kz+.61)],
F=(1-r²/R²)³. The axis is regular; all perturbation field/current traces vanish
at the radial wall. z is periodic. Density is never floored or extended outside.

`fields` analytically evaluates B, curl(B)/mu0, their r/z derivatives, the full
nonlinear E_H=J cross B/(ne), Faraday Bdot, psidot, and Poynting divergence.
Analytic density derivatives remain inside curl(E_H); using uniform n would
give a different reference. Axis limits Etheta/r and Sr/r use their exact
radial derivatives, for scalar or array coordinates. Electron pressure, flow
and resistivity are absent here.

This is an instantaneous induction/current/work verification fixture, not an
exact evolving equilibrium or formation experiment. Its finite-amplitude Hall
magnetic-energy redistribution obeys B dot Bdot/mu0+div(S_H)=0 pointwise.
The integrated poloidal and toroidal rates may exchange energy while summing
to zero. No pressure update or conservative solver state is implied.

Independent finite-difference Maxwell/derivative checks and tensor quadrature
verify the formulas before they are used to judge a numerical Hall operator.
