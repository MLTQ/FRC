# Numerical accuracy comparisons

`temporal_refinement` uses fixed spatial resolution and independent finer time
references. The smooth early-time case verifies SSP RK2. The quarter-period
case intentionally retains limiter switching and less clean ratios as a
separate measurement rather than omitting inconvenient behavior.

`reference_comparison` reruns the same wave, advection, resistive-decay and Sod
cases with the retained first-order method on identical meshes and end times.
It records both errors and full reference parameters/diagnostics so improvement
is not inferred from a historical result made with a different source tree.
