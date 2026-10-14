# Thermal finite-mass current on moving material

Extends the existing current/fluid solver with general electron/ion entropy,
their actual relative transport and its pressure work, explicitly prescribed
scalar resistance, its heat, and optional conditional Coulomb equilibration.
Two complete coupled trajectories and a time comparison are documented in
MATERIAL_THERMAL_COUPLING.md. Spatial and laboratory edge validity remain open.

Species order is electron, ion. Explicit `transport='skew'` preserves the first
runs and later high-order failure. New default `upwind` uses MaterialDualEntropy
and conservative upstream flux with matched pressure work. Its numerical
diffusion is explicit and requires spatial comparison. No heat or collision
term is removed. The full scalar-potential electric contribution
and self-drift potential -delta*mu*|w|²/(2e) are restored in `electric_material`.
`electric_field` additionally subtracts v x B at physical quadrature to give
the laboratory electric field. Induction uses the curl-free separation.
The original complete finite-mass Hall/inertia/current-stress equations remain.
Magnetic/current coefficients, fluid motion and both physical entropies advance
in the same inherited RK4 stages. No field, entropy or energy is reset afterward.

Resistance is an explicit finite nonnegative input in ohm-m, currently constant.
E_resistive=eta*J and Q=eta*|J|² use actual carrier density and drift. The heat
split defaults to mi/(mi+me) to electrons and me/(mi+me) to ions, declared
mass-weighted friction partition; callers may explicitly set the electron share.
This input is not automatically a Maxwellian Coulomb resistivity or a complete
conductivity tensor. No resistivity taper or omission of edge heat is present.

Optional `equilibration_rate` is the existing finite-mass Maxwellian Landau
temperature exchange. Qe=1.5*n*e*nu_E*(Ti-Te), Qi=-Qe. Its slow-ion/positive
Coulomb-log restrictions reject unsupported states without floors. Large drift,
neutrals, non-Maxwellian distributions and dilute edges still require regime
assessment. Turning this source off is an explicit mathematical limit, not a
default validity claim. Resistivity and this conditional source are separately
declared modeling choices, as in the existing fixed-footprint solver.

The skew option evolves local sdot=Q/(n*h), h=e*T/(gamma-1). Upwind integrates
actual heat over each dual volume and divides by its dU/ds. Both retain the
original physical heat work without clipping. Relative entropy
transport and pressure work cancel. Joule heat cancels resistive field/current
loss and exchange sums to zero. The actual mixed Faraday/Ohm rows use the same
1e-10 residual bound. Their linear refinement alters only the linear unknowns.

Geometry comes from MaterialCurrentGeometry and uses the actual entropy and
temperatures. It does not pass fake constant entropies through the isentropic
restriction. Both sides/ends still move freely with an exact empty exterior;
fields remain toroidal and flow/current meridional. Poloidal fields, moving
coil/vacuum linkage, physical heated-edge closure and FRC formation are absent.
