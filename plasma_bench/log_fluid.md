# Self-consistent motion and two physical log entropies

`LogFluidState` contains actual conservative number-density and particle-momentum
means (n vr,n vz,n r vtheta), two continuous nodal specific log entropies, and an
independent conservative cell-energy account. `initialize` sets that account to
the discrete physical energy once. Evolution never infers or resets temperature
from it. The actual kinetic metric and both thermal energies remain authoritative
physical diagnostics and their discrepancy is retained.

`derivative` supplies the single conservative continuity rate, conservative
pressure/convective stress with cylindrical hoop source, and the shared physical
bulk entropy rates. Thus momentum is now evolved, not prescribed by a callback.
The face particle integration is independently compared with the continuity
faces. `step` is atomic RK4: all stages and the final physical state must be finite
and positive without floors, hidden retries or repairs. A separate RK4 integral
of the spatial energy defect separates it from finite-step chain-rule error.

This is a smooth hydrodynamic component, not a full two-fluid plasma solver.
It has no field/current evolution, magnetic/Maxwell/drift stress, collisional
heating or Poynting flux. It has no shock admissibility theorem or diffusion.
Periodic z and impermeable regular-axis/radial boundaries are explicit fixtures.
The current mixed thermal transport cannot guarantee exact physical energy
conservation with conservative momentum; that unresolved error must converge
in independent nonlinear wave/contact tests before full moving adoption.

Measured limitation: the full-transit contact exposes a radial acoustic
instability of this baseline pressure stencil. Its source/raw failure remains
frozen and is not silently replaced by the separate radial-work experiment.
See LOG_FLUID_NUMERICS and RADIAL_WORK_NUMERICS before using this component.
