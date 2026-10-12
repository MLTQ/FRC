# Density error location and internal trace jumps

`run(density_directory, destination)` requires a complete frozen density-prefix
observation. It reconstructs all selected endpoints independently at quadrature
16/20 and requires the original reference, weight, density-change and error
fingerprints to match. No coupled solver or induction matrix is constructed.
The original summary snapshots, every state/stage hash, source/runtime closure,
derived spectra/partitions, physical fingerprints and exact time series support
a fresh complete replay. Existing source and outputs are never edited.

`calculate(case, observed)` records three complementary quantities:

- The original continuous density-error squared integral is partitioned among
  all fluid radial elements. Field quadrature cells must subdivide those
  elements, so these are exact geometric partitions of that quadrature rather
  than fractional-radius masks. The sum must reproduce the full integral.
- An FFT along the axial cell index retains each intra-cell Gauss coordinate.
  Summing all weighted mode energies must satisfy discrete Parseval. Signed
  cell modes and all spectral integrals are retained. Fractions for abs(k)<=4
  and abs(k)>=3*ncell/8 are descriptive, with no discarded modes. Cell modes
  fold harmonics separated by ncell; they are not an unaliased spectrum.
- At every internal radial interface, actual left minus right density traces
  come directly from their endpoint coefficients and are independently Fourier
  reconstructed in z. Cylindrical surface integrals use2*pi*r*dz. Surface and
  volume squared integrals have different units. Each jump norm is divided by
  the independent background density norm on that same surface, never by a
  volume norm. Initial trace jumps are retained, not removed or repaired.

Roundoff checks for the complete radial and Fourier partitions use256 machine
epsilons. These checks verify identities, not physical accuracy. Both orders
must be reported; sampled maxima need not share the integral convergence.
The original prefix is single sign with a linear continuum target, and the
different control meshes/arithmetic remain explicit. Large departures are
outside the linear target's small-amplitude regime. Trace jumps locate an
approximation feature; they do not themselves prove a failed local flux budget
or establish the mechanism of rapid growth. Nothing is filtered or fitted.

Final own-process peak memory and disk observations are separate from replay
metrics and are not continuous resource guards. Full signed trajectories,
failed-attempt replays, field/current/E and conservation audits, arithmetic/time
controls and eventual exact-vacuum formation remain separate acceptance gates.
