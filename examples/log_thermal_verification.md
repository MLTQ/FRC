# Physical thermal rate protocol

Runs all18 matched cases from `log_thermal_cases`. Both species use the same
actual conservative n/n_t and physical scalar mass, with separate gammas. Saves
log state/rate, full/reduced pressure and potential, energy projection and the
matched old lumped proxy rate. Joint physical-mass induction uses the sum of
species pressure forces and the actual density-current term.

Local thermal density/entropy and magnetic/current energy rates are independent
of the force calculation. Their combined integral equals the declared density
exchange, not zero when missing bulk work is nonzero. Exact source/reference
scope remains an instantaneous component, not full moving dynamics.

Independent quadrature12/16 reports physical pressure/temperature/proxy/rate and
matched axis norms, physical entropy and proxy defects separately. Barotropic
and isothermal controls explicitly expose absolute spurious B_t; the latter is
not asserted exactly representable. Original induction and scalar residuals
remain separate. Exceptions are preserved, make the study fail, and never trigger
silent tolerance changes. Fresh output/source archive only.
