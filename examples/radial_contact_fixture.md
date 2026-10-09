# Fixed physical radial mode on separate discretizations

Samples one supplied independent archived contact mode into the actual
conservative density, momentum and species entropies. Uses analytic poloidal
flux and cell-averaged toroidal magnetic field, with the explicit one-time
wall-conforming initial approximation. Both signed amplitudes remain actual
candidate states. No numerical eigenvector is fitted or evolved state changed.

Returns fluid, stationary base, full system, bias, fluid and magnetic tangents.
The caller retains reference identity and varies fluid/field meshes explicitly.
This opt-in fixture leaves matched_mode_fixture and all frozen studies intact.
