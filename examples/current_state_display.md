# Current state display with the latest audit status

## Purpose
Assemble a new offline 3D viewer and scientific PNG from the unchanged, replayed
trajectory display arrays. Later diagnostics have not advanced the trajectory;
the latest accepted state remains step 17 at 28.3333 microseconds. This module
updates presentation and research status without invoking a solver.

## Components
- `verify_inputs`: checks all four saved display-array hashes, the compact data,
  the completed spatial/charge checkpoint and the three complete audit replays.
- `replace_once`: fail rather than silently miss an expected frozen UI marker.
- `build`: embeds the original renderer/controller, adds camera presets and the
  current audit status, and calls `current_state_poster.render`.
- `main`: accepts a new output directory; refuses to overwrite terminal output.
  The optional server is localhost-only, with the existing bounded PNG handler.

## Contracts
The original trajectory assets, physical arrays and sources remain unchanged.
The exact embedded transformations are title, header, status, export-directory
label, and camera presets; no field values, interpolation or rendering equations
change. All source inputs and final HTML/PNG are recorded in a hash manifest.
The copied browser PNG export retains time, units and stopped-run scope. Latest
audit status is on the page and standalone scientific PNG.

The four time buttons select actual archived states. The field-line camera
presets change orientation only; the new keyboard controls support rotation
and zoom when a canvas is focused. There is no automatic time evolution.
Current density means the magnitude of the complete weak-Ampere current.
Full B includes the 5 mT bias, without exaggerated deflections.

Run end-to-end twice to new directories, compare HTML and PNG hashes, and inspect
the static image. Browser checks require an available unlocked browser. No unit
tests. Once used, preserve this source and use a new variant for later changes.
