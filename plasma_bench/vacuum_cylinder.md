# Magnetic-only periodic cylindrical exterior

The plasma occupies r<R with a fixed impermeable fluid footprint. Its outer
poloidal flux is now an evolved unknown. Full theta electric/current nodal
values include the wall; Faraday gives psi_t(R)=-R E_ind,theta(R). Ampere uses
the derivative of plasma PLUS exterior induced magnetic energy. There is no
separate auxiliary wall current, exterior density, conductivity or atmosphere.
Strong Btheta(R)=0 is retained, consistent with zero exterior toroidal field.

For the actual periodic polynomial trace N(z), form Mz=integral N N and
Sz=integral N' N'. Their discrete axial eigenmodes have k_h²=Sz/Mz. Each vacuum
mode extends as psi(r)=psi(R) r K1(k_h r)/(R K1(k_h R)). Eliminating r in [R,infinity)
gives Hvac=(2pi/mu0) Mz k_h K0(k_h R)/(R K1(k_h R)) in that eigenbasis. All axial
modes are retained. The constant mode is exactly constant flux and zero induced
vacuum field/energy. This is a Galerkin vacuum in the same polynomial axial
space, not an exact pointwise continuum Laplace solution at finite resolution.
It avoids splicing a Fourier interpolant onto a different polynomial trace.
There is no radial vacuum truncation or added numerical diffusion.

The prescribed source is an ideal periodic/infinite cylindrical azimuthal
current sheet at radius C>R, with Ktheta=B_source/mu0. Its field is uniform inside
C and zero outside. Its finite magnetic energy per axial period is retained.
It is not a finite experimental coil, real end boundary or moving free surface.
Source and induced fields superpose; vacuum induced modes extend through the
prescribed source sheet without an additional artificial wall.

`VacuumMagnetic` stores the induced field without rounding it into the source.
Unlike HarmonicMagnetic, psi_increment(R) is free. Total magnetic energy is
E_ind,plasma + E_ind,vacuum + (2pi L/mu0) B_source mean(psi_ind(R))
+ pi C² L B_source²/(2mu0). The cross term is generally nonzero. Local plasma
sampling still includes the complete B. Exterior reconstruction uses the same
axial trace and analytically decaying radial modes.

Modified-Bessel identities: https://dlmf.nist.gov/10.29 . The weak-space vacuum
elimination and source-energy formula above are derived for this module.
VacuumPlasma, VacuumInduction and VacuumGauss require these full trace maps;
the older fixed-flux geometries and archived runs remain separate controls.

Five complete coupled runs and measured spatial/time/model comparisons are in
docs/VACUUM_CYLINDER_RESULTS.md. This bounded periodic coupling is complete;
finite coils/ends, moving free surfaces and formation remain open.
