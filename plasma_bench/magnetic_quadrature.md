# Positive staggered magnetic-energy quadrature

`energy_density` integrates squared face fields with positive annular weights,
rather than squaring the interpolated cell field. Bz uses the mean of its two
axial-face squares; Br uses radial-face weights r_left/(2 r_mid) and
r_right/(2 r_mid). This equals the fluid-region part of the variational edge
energy used by the magnetic vacuum, including half dual edges on the interface.
The regular-axis face has zero weight. The approximation is consistent under
smooth refinement, not an exact continuum integral for every polynomial.

`stiffness_action` is its variational gradient divided by 2 pi/mu0. It includes
boundary reactions; its axis row is a reaction on prescribed psi=0. It does not
apply a vacuum completion or alter the field, cell total energy or pressure.
This module initially supports independent decomposition of the observed
energy mismatch. A change of evolved energy convention requires explicit
initial-state construction, pressure recovery and admissibility changes.
`selected_energy` retains the explicit historical `cell` convention and the
new `staggered` convention; invalid selectors raise. State construction,
pressure recovery, Euler admissibility and diagnostics must use the same choice.
