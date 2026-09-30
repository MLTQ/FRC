# Fluid with one shared physical radial momentum field

`FaceRadialState` stores n cell means, independent interior q=r n vr faces,
two actual cell momenta (n vz,n r vtheta), two nodal physical log entropies,
and a separately evolved conservative energy account. It cannot be loaded
from an old cell-radial state without an explicitly measured projection.

The same radial Qq field drives continuity, thermal transport, physical velocity
and kinetic energy. Axial/angular momenta retain conservative finite-volume
fluxes. Radial advection projects −Gamma dot grad(vr)+vr R n_t in the physical
variable-density radial kinetic metric. The centrifugal force is the original
physical n vtheta²/r, projected in that same metric. The metric subtracts both
species thermal face work, giving the pressure force with exact thermal pairing.

The pressure/thermal work pairing is exact. Total convective kinetic work is
not exactly conservative: axial/angular finite-volume transport and radial
strong advection use different discrete product rules. Its complete physical
energy discrepancy and the independent conservative energy account are measured.
The retained `radial_kinetic_work` method is a diagnostic of the superseded
exact radial total-work prototype, archived in face_radial_paired_initial.
Applying that work as a force gave only low-order original-equation accuracy;
it is not used by `derivative`. There is no
post-step energy/temperature reset, dissipation, floor, or numerical atmosphere.
Radial strong derivatives are cell-interior derivatives of the declared fields;
unresolved jump/shock behavior has no conservation/entropy theorem. This remains
a smooth, positive-density component requiring full original-equation testing.

All RK4 stages and final physical states validate atomically. Normal wall/axis
flux is constrained by the radial space itself; z is periodic. Magnetic/current
forces and physical vacuum/transport are not included or claimed complete.
