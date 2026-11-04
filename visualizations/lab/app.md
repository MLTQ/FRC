# app.js

## Purpose
Experiment/run navigation, live progress, launch form, actual frame replay and
diagnostic binding for the local unified app. It does not evolve physical state.

## Components and contracts
`refresh` reads `/api/state` every five seconds and adds new saved states while
preserving a scrubbed historical frame. `chooseExperiment`/`chooseRun` switch
model scope and either the IEC view or existing FRC iframe. Selection tickets
discard stale run responses. `showFrame` caches eight actual frame JSON files
and discards superseded scrubbing/playback requests.
`setup` generates unit-labelled inputs from server schema; launch/stop carry the
same-origin session token. Errors are visible and actual solver failures remain.

Playback traverses saved times only; particle count changes `Scene.count`, never
the physical model or recorded macroparticle count. Energy charts display raw
energies and uncorrected balance residual. No convergence score is invented.
Small nonzero residuals use scientific notation; they are never rounded to zero.
Run rows show actual duration and integration timestep to distinguish controls.
Selected run inputs are displayed alongside model scope. The absolute residual
is also normalized to current particle kinetic plus cumulative injected/absorbed
kinetic energy, so large electrode charging energy cannot hide particle error.
Selecting a preparing run clears the previous run's canvas, timeline and metrics;
its data are never shown under the new run's status.
Run links retain `?run=ID`; an existing FRC iframe's run selector updates the
outer library selection. Rows show source conditions as well as duration/dt.

HOPS runs retain their explicit mode/label and show actual microseconds,
live-ion RMS radius/speed and coherent/residual motion energies. The additional
radius/voltage chart uses saved state only. Saved inputs hide legacy sinusoid
controls on feedback runs and expose their actual sensor/response thresholds.
Externally coordinated HOPS/control runs have no misleading app-owned Stop
button; their scientific campaign retains its own STOP mechanism.
HOPS is expanded as Harmonically Oscillating Plasma Sphere in the heading;
small nonzero diagnostic values use scientific notation rather than display zero.

`FieldsView` adds read-only native E/current/density/phase views for every saved
IEC frame, with an optional explicitly partial B reconstruction. Playback waits
for both raw and derived data, then commits all panels together in one browser
animation frame. `frameLoading` prevents overlapping playback requests; old
frames remain visible while scrubbing loads. Run changes clear old
fields and invalidate stale requests. The first metrics prioritize ion retention
and core population; electron populations remain explanatory field diagnostics.
