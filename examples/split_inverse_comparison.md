# Opt-in arithmetic baseline and amplitude controls

Requires the same-grid old and new actual inputs and complete signed fluid
rates to match bitwise. Magnetic/current/electric results may change through
strict arithmetic refinement; saves their centered coefficient differences
and reports Euclidean norms explicitly distinct from physical error norms.
The full physical continuum error norms remain in the rate archives.

For doubled radial field resolution, compares amplitudes 1e-4 and 1e-5 on the
same grid. Checks every original strict acceptance bound, input hash, terminal
status and quadrature16/20 rate/initial volume/axis/wall norm difference, with
absolute scales retained for tiny denominators. Does not assume all meaningful
digits are physically accurate. Fresh full replay required; no failed input
is accepted here and no full trajectory or formation gate is closed.
