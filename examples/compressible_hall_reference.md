# Compressible Hall continuum reference

Independent cold-electron, massless-electron Hall waves at finite ion sound
speed. No production Hall spatial/time operator is imported. The six amplitude
equations follow continuum momentum, ion adiabatic pressure and Faraday laws;
equal-energy coordinates make their generator anti-Hermitian. `eigh(i A)` gives
signed real frequencies and polarizations, not an assumed incompressible curl
helicity. `dispersion_terms` independently eliminates the amplitudes to the
sixth-degree dispersion relation. A q=0 generator is a Cartesian parallel
limit; the cylindrical class always uses a positive J1 root.

`fields` returns perturbations (including current and linear electric field),
with SI vectors in r/theta/z order. Density is p1/cs² for these isentropic waves;
the independent zero-frequency entropy mode is not included. `amplitude` is
the norm of equal-energy coordinates in m/s. Branch 0/1/2 means ascending
positive frequency for the chosen coefficients, not a branch identity through
crossings. Signed bias, axial mode and frequency are all retained. `hall=False`
is an explicit MHD limit, not a zero carrier density. Invalid inputs fail.

Pressure belongs entirely to adiabatic ions (pe=0). This avoids silently
omitting electron enthalpy transport from a warm barotropic closure. Fully
ionized positive uniform density, no resistance/electron inertia, conducting
slip radial wall, periodic z. Resistive Jz at the wall is generally nonzero:
these fields do not satisfy the scalar-resistive conducting boundary. Linear
fields are not exact finite-amplitude nonlinear solutions. Energy is the
quadratic perturbation invariant, not the large background total energy.

Tests must check physical Maxwell/species derivatives, radial traces, volume
energy, the independent dispersion and limits; matching the matrix alone is
insufficient. Derivation and forthcoming nonlinear integration gates are in
COMPRESSIBLE_HALL_REFERENCE and NONLINEAR_HALL_PLAN.
