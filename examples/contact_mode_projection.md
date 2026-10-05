# Matched Fourier/Bessel physical projection

Projects the independently reconstructed production perturbation onto physical
J0/J1 radial modes and paired Fourier harmonics. Radial inner products use the
same cylindrical measure as the fields; analytic Bessel norms are pi R² J0(qR)².
The radially constant scalar/axial mode is included. Even-grid Nyquist is omitted
to retain conjugate pairs and a real band-limited perturbation; its lost norm is
included in the captured-versus-original norm diagnostic. This is analysis only,
not a filter or modified production state.

The frozen independent continuum parent operator acts on every projected mode,
including density, all velocity components, both specific entropies and B.
Its radially constant reduction follows continuity, total-pressure momentum
and both entropy advections, with no transverse J1 mode. All parameters are the
specified argon contact fixture, not general laboratory defaults. Report each
component's captured norm and eigen-equation residual, varying Fourier/radial
cutoffs and physical quadrature. Uncaptured norm, base interpolation error and
boundary/face mismatches must remain explicit in interpretation.
