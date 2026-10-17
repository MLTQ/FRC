# Ground-molecule total ionization

`MolecularBEB.from_nist` loads attributed H2/N2/O2 orbital constants and verifies
their source hashes. The orbitals are B (binding eV), U (kinetic eV), N (occupied
electrons). This implementation deliberately accepts only the supplied Q=1,
unmodified n=1 orbitals, without a double-ionization multiplier. D2 is unsupported.

For t=E/B, u=U/B and S=4*pi*a0²*N*(R/B)², each orbital contributes

    sigma = S/(t+u+1) * [log(t)/2*(1-1/t²) + 1-1/t-log(t)/(t+1)]

above threshold and exactly zero at/below threshold. The code uses log1p and
algebraically equivalent differences near threshold. NIST CODATA 2022 values
a0=5.29177210544e-11 m and R=13.605693122990 eV are explicit. Sources:
[Bohr radius](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0),
[Rydberg energy](https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev).
These replace the rounded explanatory constants in the first comparison;
they are independently sourced values, not fitted amplitudes. Output is SI area in m².
The supplied orbital contributions sum to the declared total; no new channels
or product probabilities are inferred. The explicit energy range is 0–5000 eV,
matching the independent downloaded tables; values above it reject.
Scalar and arbitrary array energies preserve their input shape; orbital output
adds one trailing dimension and total cross section sums that dimension.

`maxwellian_integral` assumes stationary targets and isotropic, zero-drift
Maxwellian electrons. It returns the integral only within that energy range:

    k = sqrt(8*e/(pi*me)) / Te^(3/2) * integral sigma(E)*E*exp(-E/Te) dE.

Threshold-relative coordinates and geometric subintervals resolve the tail with
Gauss–Legendre quadrature. The probability above the supported energy range is
reported separately; the electron distribution is not renormalized or silently
extended. That probability is not an upper bound on the unknown rate outside
the range. The result key explicitly says `rate_within_energy_range_m3_s`.

This is a cross-section/rate component, not a molecular discharge source.
Total ionization does not determine intact/fragment products, momentum sharing,
internal excitation or energy reservoirs. It therefore cannot replace the argon
reaction network or turn the old H2 mass/inventory switch into molecular physics.
The NIST tables independently verify numerical implementation of their theory;
they do not constitute an independent experimental validation of the theory.
