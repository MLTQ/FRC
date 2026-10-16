# Actual pulse versus bias-only response

Reads two completed configured experiments with matching tube, initial plasma,
mesh and duration. Plots the prescribed centre field, actual boundary radius,
azimuthal current and accumulated Joule heat. A dashed radius curve follows the
material midplane, interpolated at reference z=0 on the actual side nodes. This
distinguishes central compression from outer expansion; it is not a density or
confinement proxy. The configured symmetric cases retain their physical midplane.
Final signed Bz and poloidal-flux
sections are reconstructed from the saved finite-element fields at common
material labels, each mapped to its actual physical coordinates. Both field
panels use one signed color range from their recorded field extrema, rounded
outward to two significant figures; the radial visual scale is explicitly
expanded. No state is evolved or adjusted and no interpolated cloud substitutes
for the finite-element fields.

Purpose: distinguish the pulse's induced magnetic response from the large free
expansion already present in the control. Reversal is not confinement or an FRC
formation certificate. Cold-edge/transport, initial ionization assumptions and
local current resolution remain limitations of the input simulation. Output is
a scientific PNG/PDF and a concise comparison JSON with source/summary hashes.
Pulse amplitude, time and chamber/gas labels come from the saved configuration.

Also compares two pulsed coil arrangements: labels identify pulse strength and
coil count from each input, rather than assuming the first run is bias-only.
Particle-weighted RMS radius histories use the actual Q1 geometry and fixed
particle inventory; dotted curves distinguish bulk extent from surface radii.
The heat panel shows Joule heating and explicitly numerical compression heating
separately in mJ. Missing numerical heat in older archives remains null and is
not drawn. The JSON retains both heat histories, declared viscosity coefficients,
minimum radii and maximum observed external-work residual relative to initial
energy. Neither cumulative heat integral repairs the trajectory's energy budget.
At peak drive and final time, drift/thermal thresholds 0.1 and 1 are reported as
fractions of volume, particle inventory, electron thermal energy and integral J²dV.
The last weight is proportional to Joule power for these uniform prescribed-eta
runs. These content measures distinguish a dilute-edge maximum from bulk impact;
they are not a declaration that other samples have valid transport.
At each run's recorded peak drive, each group's independent
exterior flux linkage is compared with its volume-current linkage. The finite-
space discrepancy is reported unchanged; it never repairs work or currents.

The same tool also compares different physical collision closures on otherwise
matching plasma, mesh, timestep, integrator and numerical compression treatment.
Coil groups and/or collision models may change; every change is named in JSON
and the legends. Prescribed eta versus paired Maxwellian heat is labeled as
collisional heat, with numerical ion heat separate. Stored signed species heat
is retained when available, otherwise null. Peak/final response now includes the
actual collision force, so drift-content weights also use J.E_collision; J²
alone is not a heat weight for nonlinear finite-drift collisions. Calculations
are saved before plotting. Comparison of two closures is not a validation of
either closure, and two simultaneous physical changes cannot be attributed to one.
