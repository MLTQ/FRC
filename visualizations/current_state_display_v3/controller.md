# Signed current controls and export

Extends the original linked camera, four saved times, per-frame scaling, optical
opacity and PNG export with `|J|`, `Jr`, `Jtheta`, `Jz` selection. Every selection
updates the displayed quantity, current readout, symmetric signed legend, ramp
and explanatory text together. Shared limits span all four actual snapshots.
No time interpolation is performed. Camera presets remain a separate frozen
extension. Export labels and gradients follow the selected component and frame.

Consumes the builder's component ranges and four-channel current textures.
Signed color conveys the cylindrical component sign, not a Cartesian vector
direction or an inferred current loop. Export retains the failed-duration and
verification-fixture status. No physics or saved states are changed. Preserve
once used; end-to-end browser verification is separate from PNG rendering.
