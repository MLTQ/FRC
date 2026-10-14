# Complete heated material expansion

Runs the actual MaterialThermalCurrent equations from the same compact physical
initial density/temperature and analytic B/J profile as material_current_expansion.
Positions, momenta, independent field coefficients and both species entropies
advance together for one complete acoustic crossing R/sqrt(10e/(mi+me)).
No motion, temperature or field shape is prescribed after initialization.

The selected constant scalar resistivity is a declared model input. Joule heat
uses the actual J² and remains present at all material quadrature, with the
declared mass-weighted electron/ion share. The existing conditional Maxwellian
temperature exchange also acts on the evolving states; no Coulomb-log floor.
This does not claim that constant resistance plus this exchange is a complete
laboratory transport closure or that the cold heated edge is physically valid.
No exterior material, poloidal B, external coils or ionization is present.

The independent scalar virial remains I''=4E-2Wb for gamma=5/3, now with changing
species entropies and actual dissipation converted to internal energy. The
driver integrates the magnetic term separately without feeding discrepancies
back into the solution. Initial velocity is zero. Joule/exchange integrals use
trapezoidal time quadrature and are observational, not exact RK4 stage integrals.
All accepted canonical states, actual species energies/temperatures/entropies,
physical equation residuals and source closure are saved. No unit tests are run.

The existing isentropic expansion is left intact while its higher-order run is
live. This driver exercises actual thermal coupling rather than another isolated
operator check. Both use the same analytic compact initial field/current target.

`--transport skew` retains the original spatial operator and its preserved
high-order failure. New default `upwind` uses explicit dual-volume diffusion
with matched pressure work. Minimum species entropy and actual outflow rate
are saved; no entropy values are clipped to their initial bound. Uniform initial
entropy is representable in either method. No heat/exchange/threshold is changed.
