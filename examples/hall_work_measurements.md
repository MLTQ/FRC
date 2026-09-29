# Hall spatial measurement records

Runs every declared case and preserves source/configuration identity, complete
corner/cell initial fields and rates, independent reference rates, projected
electric coefficients, carrier quadrature and local Poynting flux/rates. Each
lossless NPZ has a verified hash. These are instantaneous spatial audits; no
time history, accepted evolution or formation result is implied by completed.

Physical-volume quadrature compares actual electric-space current and actual
face-interpolated magnetic rate with continuum fields. Conservative g means,
physical r²-weighted toroidal rates, axis rates, corner psi rates and local
energy defects are separate. Zero-reference rates have null relative error and an
explicit absolute local-energy defect; they are not divided by an arbitrary
floor (the initially zero-toroidal case has zero instantaneous magnetic work).
Raw perpendicular electric samples are distinct from their projection. Global
cancellation is normalized by |J| |E| work.
Reference volume quadrature is checked at orders 10/14. Common electric-space
quadrature and carrier sampling have separate protocol controls.

Opposite perturbations about the same uniform bias algebraically isolate the
quadratic Hall term, but their floating-point subtraction becomes ill-conditioned
at small amplitude. A direct Hall(b,b) evaluation on the pure perturbation
separately resolves this term without subtracting O(bias*perturbation) rates.
Both results and their discrepancy are retained. This is used for amplitude
scaling, not as an exact nonlinear time reference. Linear Hall rates retain
finite-amplitude residuals. No failed
axis accuracy gate is overwritten by successful conservation measurements.
