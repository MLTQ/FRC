# Moving plasma field display shell

Offline visualization of actual saved full-vector material states. The title,
scope notice and source description distinguish initialized closed flux from
formation and distinguish axisymmetry from 3-D stability. Two linked canvases
show magnetic paths and sampled current density. No hardware controls or
simulation parameters are hidden behind display controls.

DATA and VIEWER markers are filled by material_vector_display. Every slider
position is one actual state; playback does not interpolate time. Current magnitude
is the default, with all signed components available. Color-scale
mode is explicit and field/current values retain physical units. The moving
outline is the material boundary, not a proposed experimental chamber.

The timing label distinguishes the numerical integration step and total steps
from the selected display states. Frame gaps come from actual saved timestamps.
Previous/next buttons inspect one frame at a time; the1/5/10/20fps selector changes
wall-clock playback speed, never simulated timestamps or solver integration.
Dense local displays load their actual saved frames on demand through the viewer.
