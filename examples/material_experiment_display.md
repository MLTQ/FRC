# Actual biased-coil experiment display

Reads a completed configured experiment through its unchanged physical builder.
Eight evenly spaced saved states plus unipolar transition times show the initial
bias, pulse, and late relaxation. The view never interpolates in time or evolves
new states. Every selected raw array and source summary has an input identifier.

Field-line paths use actual plasma B and the actual deformed vacuum B, plus exact
prescribed coil fields. Linear physical-space interpolation and finite tracing
steps in both directions are display approximations. Axis/interface errors and
near-filament interpolation are not repaired or turned into a topology claim.
The field legend covers the traced paths and plasma samples, not the divergent
ideal filament. Current clouds use actual en*w, including all components, only
inside plasma. No vacuum current or exterior matter is drawn.

Grey tube outlines, copper coil loops and moving plasma surfaces use one physical
length scale. The scientific still selects the largest-current saved state;
the offline HTML exposes every selected time, linked rotation/zoom, signed current
components and explicit color-scale controls. Source field at the vacuum centre
is labeled separately from the actual plasma field. The long tube has a spatially
nonuniform field even with Helmholtz-spaced coils.

The visible scope preserves assumed preparation/ionization, prescribed transport,
coarse mesh/current limits, and the distinction between field lines, topology and
particle confinement. This is a visualization of the simulation, not a new
validation or reporting framework. The older compact-flux display is preserved.

A second scientific figure samples the actual material basis at three saved
times. Signed Bz and Jtheta use physical moved coordinates; flux contours show
the actual total psi, with zero highlighted only if crossed. Dense plot sampling
adds no simulation accuracy and is not a separatrix or closed-flux certificate.
Unlike the 3-D view, the 2-D sections explicitly expand the radial scale for
readability; both axes retain metres and the figure states that convention.
The initial camera zoom fits the long tube while preserving all proportions.
