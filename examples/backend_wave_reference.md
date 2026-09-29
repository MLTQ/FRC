# External-engine continuum references

`ParallelWave` derives the positive frequency of both cold circular Hall branches
from the massless-electron two-species quadratic. B=(A cos phase,-h A sin phase,B0)
and transverse ion velocity is -k B0 B_perp/(mu0 n mi omega). Constant density,
zero thermal pressure, no collisions, periodic z and physical proton mass apply.
This is the common cold transverse limit of hybrid PIC and our linear Hall
component, not an argon formation experiment. k di=1 exercises Hall dispersion.
The magnetic amplitude is 1e-4 of the background unless explicitly overridden.

`vacuum_field` is an independent source-free Maxwell plane wave with Ex=c By.
No plasma or density floor is part of this reference. It verifies the full
electromagnetic solver, not hybrid vacuum or plasma-edge coupling.

`relative_error` and `mode_amplitude` measure physical samples on their actual
staggered positions. They do not import WarpX or our production operators.
Consumers: `warpx_comparison` and `backend_comparison_report`.
