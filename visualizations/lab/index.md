# index.html

## Purpose
Unified local experiment UI shell, with experiment/run navigation, an honest IEC
workspace, existing FRC viewer, and physically labelled new-run inputs.

## Contracts
`app.js` populates registry, manifests, frames and diagnostics. `scene.js` draws
actual saved 3D points/wires and a physical potential cut. Display-particle and
time controls never alter the solver. FRC is an existing-viewer iframe adapter.
Model scope remains visible alongside the visualization and in the launch form.
No third-party scripts, generated scientific data, or external analytics.

Saved-input details show the actual selected run configuration, including its
source/drive parameters; opening a new-run form does not substitute defaults.
HOPS and its controls add a saved radius/voltage chart with explicit separate
scales. The panel is hidden for legacy IEC and FRC views.

Additional `fields.js` panels show 3D E/current/optional plasma-only B, native
particle moment slices, full ion radial phase space, radial density averages and
live population history. Their scope and physical units stay visible. B is not
an evolved total field. Ion density is the default scalar, with ion retention
and core count prioritized in the metrics.
The module entry and its imports use the same explicit asset revision so a
normal reload adopts a coordinated viewer update rather than cached old code.
