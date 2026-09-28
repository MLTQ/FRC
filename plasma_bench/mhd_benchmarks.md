# Independent dynamic verification fixtures
`radial_fast_mode` uses analytic Bessel eigenfunctions and the ideal-MHD fast
speed, with volume-averaged density/pressure/velocity. It is a linear solution
accurate to O(epsilon²), not derived from the numerical update. A radial J1 zero
satisfies the slip conducting wall. `resistive_bessel_mode` supplies a weak-field
diffusive eigenmode whose fluid feedback is asymptotically negligible.
`advected_flux` supplies a weak 2-D flux loop advected along periodic z.
`sod_initial/sod_exact` provide a gas shock tube and its independent rarefaction,
contact and shock solution. Compare away from periodic boundary interactions.
These fixtures test algorithms; none represents a demonstrated formation event.
`radial_pressure_balance` supplies an axial Bz=B0(1-alpha r²), whose current and
radial magnetic force are analytic. Optional p+B²/(2mu0)=constant gives a known
static equilibrium. Evolving it measures spurious flow from numerical diffusion;
it is a straight-field equilibrium, not an FRC equilibrium or formation result.
`longitudinal_rarefaction` separates two gas streams parallel to a uniform axial
field. The field exerts no force and must remain unchanged. Its independent
`double_rarefaction_exact` uses Riemann invariants for two isentropic fans and
a low-density/low-pressure central state. It explicitly excludes initial data
that produce true vacuum. Compare only before periodic-boundary shocks arrive.
