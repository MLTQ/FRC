# Independent density observations of accepted prefixes

`run(manifest_path, destination)` reconstructs every saved density in each
declared prefix, at independent quadrature 16/20. The immutable input manifest
contains the original summary snapshot, directory and selected prefix length.
The live summary may advance; all selected initial and step bytes must remain
identical. The source closure, independent continuum input, exact raw hashes,
physical fingerprints and all measured time series are retained for fresh replay.

This lightweight observation builds the fluid node layout and polynomial/Fourier
sampler only. It constructs no coupled plasma solver, induction matrices or time
integrator. Final own-process peak memory and disk observations are stored
separately from scientific replay metrics; these are not continuous guards.

`calculate(case)` samples the actual scalar density by independent root-product
polynomials and Fourier sums. Cylindrical Gauss weights integrate over field
cells, which must subdivide the fluid elements. The continuum target is
`n_base(z) + epsilon * exp(lambda*t) * dn_reference(r,z)` using the unchanged
radial1/N64 analytic Bessel/Fourier reference. The independent density error is

    sample(n(t)-n(0)) - epsilon*expm1(lambda*t)*dn_reference
      + [sample(n(0))-n_base-epsilon*dn_reference].

Subtracting saved coefficients first avoids unnecessary background cancellation;
the initial representation error is explicitly retained. Each time reports the
error relative to the growing linear signal, actual density change over the
initial signal, nodal and sampled minimum density, and maximum fractional
change relative to the initial density. No unlike-unit fields are combined.

This is a single-sign density observation of an accepted prefix. It does not
require or assert terminal trajectory success, measure all physical fields or
conservation, remove finite-amplitude effects, fit a growth rate, identify an
unstable eigenmode, or establish the source of any growth. A linear target is
not an exact nonlinear solution at large departures. Prefix-observation replay
does not replace full trajectory state/stage replay. Existing full signed
audits, time/amplitude/mesh controls and eventual exact-vacuum formation gates
remain required. No source, simulation state, tolerance or diffusion is changed.
