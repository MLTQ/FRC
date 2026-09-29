# Independent electron reference fields

`ElectronReference` differentiates smooth cylindrical density, temperature and
toroidal magnetic fields analytically. It supplies independent Biermann gdot,
electron entropy drift, current and physical pressure electric field, including
the finite axis limit. Density and temperature perturbations have different
radial profiles and phases so their gradients are not parallel. Wall density
is constant, toroidal field vanishes, and battery rate vanishes at the wall.

`cell_average` independently integrates exact fields in physical annular
volume with tensor Gauss quadrature. It does not use discrete current, energy
or pressure-curl routines. `drift_pressure` is a true rigid electron translation:
uniform n and Btheta=slope*r give constant Jz and ve,z=-2*slope/(mu0 n e).
Electron pressure is transported although the bulk velocity is zero; pressure
Ohm is curl-free and the toroidal field stays fixed. Nonzero constant wall Btheta
must be supplied, rather than inventing a surface return current.
