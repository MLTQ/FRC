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
particle inventory. At each run's recorded peak drive, each group's independent
exterior flux linkage is compared with its volume-current linkage. The finite-
space discrepancy is reported unchanged; it never repairs work or currents.
