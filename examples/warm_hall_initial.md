# Warm Hall initial conditions

Linear-wave fluid variables are independently continuum-averaged by the existing
wave reference initializer, with exact node psi and pe0+pe1. No coarse-cell
energy adjustment makes recovered species pressures artificially exact.

Nonlinear cases use the configured analytic poloidal flux, toroidal g, positive
cell density and node pe, with explicit wall-varying density/pressure factors.
Total pressure is the Q1 electron average plus specified uniform ion pressure;
initial compatible energy is constructed once from that explicit state. The
case starts at rest and is unbalanced, not an exact equilibrium/formation
solution. Subsequent total energy changes only through fluxes. All coefficients
are in the case configuration; geometry, ion mass, gamma and tolerances are
passed explicitly to the coupled solver.
Both iteration budgets and all fluid reconstruction/CFL choices are passed
from the saved case instead of relying on constructor defaults.
