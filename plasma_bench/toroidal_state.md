# Compressible state with toroidal degrees of freedom

`ToroidalState` has a separate explicit six-variable layout rho, mr, mz, total E,
axial angular momentum density L, and g=Btheta/r, plus corner psi. It is not an
implicit reinterpretation of `MHDState`. `primitive` returns rho, vr, vz, p, Br,
Bz, Omega=vtheta/r and g. SI units; all pressure/density failures raise.

`radius_squared` is the exact annular volume mean of r². Recovery assumes a
cellwise regular vtheta=r Omega, Btheta=r g, constant rho/Omega/g. Thus
L=rho <r²> Omega and the toroidal kinetic/magnetic energy use <r²>, including
the axis cell. This is exact for uniform-density rigid rotation and Btheta=c r;
subcell variations introduce a declared second-order modeling error. It does
not pretend to recover arbitrary subcell energy from averaged momentum.

`from_primitive` constructs this numerical energy convention; exact continuum
fixtures may instead supply independently integrated conservative variables.
Diagnostics retain total angular momentum, physical toroidal flux (integral
Btheta dr dz), all energy components and CT divergence. No energy reset, thermal
floor, Hall or electron-current energy is included in this massless-MHD state.
