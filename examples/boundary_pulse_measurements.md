# Driven boundary comparison

Runs the unchanged two-coil fixed-chamber pulse with the compatible RK3 method,
with or without the explicit first-order wall/end strip. The comparison changes
only that declared selection. Both use the same physical resistivity, initial
gas, prescribed sources, numerical interface resistance and energy convention.
No formation or physical argon interpretation is attached to these unit fixtures.

`protocol` fixes grids 16/32/64/128, half CFL and source time cap at 64, and a
farther magnetic matching surface at 32 for both methods. `measure` records
every output time's budgets and saves full fluid conserved/primitive fields,
CT total flux, stiffness current and numerical/physical interface eta ratio at
0.04, 0.06 and 0.10 s. The NPZ is lossless and has its own SHA256 in the JSON.
The field shapes and units are explicit; current and eta live on fluid nodes.
Primitive fields retain exact samples rather than interpolated plotting data.

Engine, common measurement, probe, driver and dependency versions are recorded.
`run` reuses only completed matching records and verifies the associated NPZ.
Progress is separately marked incomplete, and failures retain accepted state.
Comparing source work alone is insufficient: the associated report must measure
bulk and boundary profiles, kinetic response, mesh/time/matching sensitivity and
physical versus numerical heating. Adaptive limiter counts exclude the fixed
strip, whose support is separately recorded.
