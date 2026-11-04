# fields.js

Read-only field and phase-space panels for the full native IEC/HOPS state.

`FieldsView` requests `/api/observables`, caches eight responses, and serializes
requests while coalescing scrubbed frames. Selection generations discard stale
responses. Loading preserves the last complete view. `prepare` resolves a
native state without painting; `commit` updates derived panels together with
the particle frame inside one animation-frame callback. Queued scrub requests
supersede old requests, which resolve without repainting. Playback waits for the field
request so lines and particle panels refer to the same saved time. Mode changes
retain the previous labelled view until their replacement is ready. The display
particle-count slider does not affect these full-population observables.
Field-mode changes go through the app's frame coordinator using the requested
timeline index. They cannot revert an in-flight scrub to the old displayed time.

Electric lines, normalized current glyphs and optional plasma-only magnetic
lines share an interactive `Scene`. Magnetic output is prominently labelled as
a regularized partial diagnostic, with missing currents/steel/force physics.
It is never presented as the evolved total magnetic field.

The y=0 selector shows total/electron/ion |J|, ne, ni, signed rho or |E|. Its
nonnegative colors use log mapping with a labelled sampled maximum; charge uses
a symmetric signed linear scale. No spatial smoothing or invented frames.

Ion phase space uses all saved ions; spherical radial profiles are explicitly
diagnostic averages. Live-population histories use log10(1+count), preserving
zero and small populations. Every plot has physical units and a time marker.
Frame captions retain up to six decimal places in microseconds, so 1 ns
archives are not rounded into indistinguishable 0.01 µs labels.
