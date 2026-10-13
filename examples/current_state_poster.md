# Latest scientific field poster

## Purpose
Render the unchanged final accepted field arrays with the latest research status.
This is a separate presentation variant of the frozen trajectory poster; its
data and geometry remain the same, with a closer camera and updated captions.

## Components
- `render(source, summary, data, target)`: plots the supplied last frame and saves
  a 2880 by 1872 PNG. `current_state_display.build` verifies its inputs first.
  Current magnitude is evaluated from all three archived vector components.
  Magnetic paths retain every archived point, full bias and actual deflection.

## Contracts
No solver, resampling of physical arrays, state modification or new time sample.
The common linear scales cover all four saved times. Surface subsampling,
alpha, camera zoom and line spacing are presentation choices. The coordinate
scales preserve the physical cylinder aspect ratio; axis labels use millimeters.
Current is shown as cutaway surfaces, rather than the interactive optical volume.
Magnetic lines use total |B| in mT, and current uses kA/m².

The state is an axisymmetric filled periodic verification fixture, not an FRC.
The image states both the step-18 rejection and unresolved rapid growth.
New audits do not advance the last accepted time of 28.3333 microseconds.
Execute complete rendering, inspect the PNG, and repeat to a new directory for
byte comparison. No unit tests. Freeze this source after use.
