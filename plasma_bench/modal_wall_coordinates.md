# Separable admissible physical coordinates

Uses the same fixed mean density, sound speed and equal-energy magnetic scales
as PlasmaCoordinates, removing the same inactive axis and normal-wall fluid
slots. A small radial trace-null basis is tensored implicitly with axial identity.
There is no dense full-grid null-space matrix. Both psi and toroidal stages are
constructed admissibly; no accepted state or derivative is repaired.

Coordinate ordering is radial/component row followed by axial point. `pack` and
`unpack` remain real physical linear coordinates. `fourier`/`physical` transform
complex *preconditioner* vectors with an unnormalized FFT/inverse FFT. Physical
nonlinear residuals still use the real coordinate infinity norm; the radial
orthonormal basis and its ordering are explicitly part of this new norm contract.

Matching axial field/fluid counts are required for the exact uniform modal
inverse. Other dimensions remain independent. Actual nonlinear fields, density,
current and both species are never Fourier filtered or replaced by a background.
