# Moving nonuniform Coulomb temperature contact

Advance the complete CoulombExchangeWallPlasma for one axial transit100 us:
R=0.035 m,L=0.24 m,n=1e20 m^-3,Ar+,both gamma5/3,vz=2400 m/s,Bz=5 mT.
Initially Te=5+0.3 cos(2 pi z/L) eV and Ti=6-Te. Density and total pressure
are constant. The exact reduced-model solution has constant n,v,B and zero J;
Te/Ti relax along z-vt characteristics. The axial electrostatic pressure field
is solved by the original species/current equations, not externally forced.

The reference calls the separately integrated scalar temperature clock from
coulomb_temperature_relaxation for each characteristic's initial temperature.
No production collision/RHS routine generates a later reference state. The
nonuniform state advects and thermalizes inside every original Gauss stage;
no operator splitting, prescribed evolved profile, damping, field filter,
energy reset or temperature floor is used.

F4,B8 radial spaces stay fixed because the exact solution is radially uniform.
Run axial16/32 with the same80 steps to assess actual axial resolution effects.
The method uses its original1e-12 absolute/1e-10 relative acceptance bounds,
sixth-order centered axial fluid transport and degree5 compatible field basis.
Saved initial errors expose temperature interpolation. Each accepted state and
stage is retained, and independent q20 reconstruction measures temperatures,
density/pressure, velocity/current/B departures, energy and entropy. Conservation
changes reference the actual initial discrete state. Native and physical energy
remain distinct; density error is against the analytic n, while particle change
uses coefficient subtraction before reconstruction.

Physical entropy is sum(C_s)/(gamma-1) in kB units; its exact integral changes
through local thermal equilibration. Report entropy change and error separately.
All species share the same bulk velocity analytically; common advection is not
relative electron-ion drift. The measured current/thermal-speed ratio exposes
departures from that assumption.

This idealized component benchmark omits heat conduction, physical friction,
neutrals, ionization and radiation. Its nonuniform temperature profile should
not be interpreted as a complete collisional laboratory prediction. Boundary
conditions remain periodic z, impermeable wall, fixed psi and zero Btheta.
It establishes no exact-vacuum interface, formation or resolution of frc-7mm.7.

Resources are observed after setup and each accepted step (own4 GiB RSS,
256 MiB free disk); this is not a continuous native-call guard. Existing source
and result helpers preserve every endpoint. No unit tests, new verification
framework or automatic replay is added.

All three axial16/80-step, axial32/80-step and axial32/160-step transits
completed. Axial refinement reduces final Te/Ti errors by50/55.7 times;
time refinement reduces energy drift by15.99 times. Fine temperature errors
remain below1e-6 eV. Tiny current/B departures do not converge with time
refinement. Read docs/ADVECTED_COULOMB_RESULTS.md for evidence and limits.
