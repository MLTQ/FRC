# Independent electron pressure energy oracle

Extends the analytic `ElectronReference` with physical Ep cross B/mu0 and
-gamma/(gamma-1) (pe/n) J/e, plus local d(epsilon_e+HB)/dt from independently
differentiated thermodynamic and Faraday equations. This module imports no
production current/energy transport. Gauss integration gives physical primal
face powers and cell-average energy rates in SI.

The optional adiabatic barotrope pe proportional to n^gamma has constant
specific entropy and zero pressure induction; total local energy rate is
exactly zero. Its generally nonzero J dot Ep supplies a physical work scale
for an absolute cancellation error; a zero-reference relative norm is not
invented. Uniform-density and nondefault-gamma controls use the same formulas.
