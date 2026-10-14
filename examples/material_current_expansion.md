# Moving finite-mass current trajectories

Two concrete end-to-end runs use MaterialCarrierCurrent. `reference` prescribes
homologous meridional expansion lambda=1+t/tau, uniform initial n, constant
species entropies and an initially weak Btheta=B0 J1(kr r) cos(kz z), with
J1(kr R)=0 and kz=pi/(2L). The physical electron/ion masses set
chi=mu*(kr²+kz²)/(mu0*e²*n0)=1. Linear extended MHD independently gives
b(t)/b0=(1+chi)/(1+chi*lambda), Bphysical=b/lambda², Jphysical=J0*b/lambda³.
Current kinetic energy is comparable to magnetic energy; ideal flux freezing
would give a different result. Full nonlinear Hall/self-drift terms remain in
this weak-amplitude reference. Motion is externally prescribed only in this
explicit verification mode, and field/current work is reported as such.

`free` evolves the actual simultaneous material positions/momenta and independent
magnetic coefficients. n0=nc*[(1-r²/R²)(1-z²/L²)]^(3/2), T0=Tc*n0^(2/3)/nc^(2/3)
sets equal constant species entropies; pressure/density vanish on side/ends.
The initially compact toroidal field is projected into the actual boundary
space. There is no external force, prescribed shape, coil or exterior matter.
Hall, finite current inertia, actual geometry terms and current stress remain.
The full source-free XMHD virial identity is I''=4E-2Wb for gamma=5/3,
including current kinetic energy in E. No heat or energy residual is repaired.

Both modes retain every accepted state and source, original direct-solve residuals,
physical energies/current/drift and resource guards. The constant-entropy,
zero-resistance invariant sector is a mathematical model comparison, not a
laboratory edge/transport validation or FRC formation. General entropy/current
transport, heating, poloidal fields, vacuum field and finite coils remain needed.

Initial projection is now explicit: `magnetic` preserves the earlier input
approximation; `magnetic_current` (new default) minimizes physical B plus carrier
velocity error against the same analytic initial fields. This addresses the
measured85%/29% initial current-energy excess from magnetic-only projection.
No continuing state, heat or energy is corrected. Archives identify their source
and new records name the projection. Independent initial B and J errors remain
in each record; improved total energy alone does not validate the representation.
