# Static reference figure and independent integrals

Renders analytically prescribed Bz, pressure and toroidal current with closed
flux contours, spherical separatrix and magnetic null. Every panel is labeled
as an equilibrium reference, not evolution or formation. Parameters are unit
verification constants. The exterior pressure and current are exactly zero.

Independent spherical Gauss quadrature checks continuum pressure/interior field
and full-space induced energies, using x=a/s for the infinite exterior. Saves
analytic targets, measured quadrature errors, full reference parameters and
separate source fingerprints. Output PNG requires visual review. Density,
temperature, physical transport, lifetime and stability are not assigned.
The compact equal-aspect layout keeps spherical geometry undistorted; the
current palette makes zero-current vacuum dark while retaining signed values.
The existing topology detector is independently checked against the exact
O-ring, two axis nulls and spherical separatrix volume on four refined nodal
grids. Complete detector output and volume error are saved, with detector,
geometry/operator hashes and library versions. Sampled largest closed contours
are distinguished from the explicitly axis-attached separatrix estimate.
