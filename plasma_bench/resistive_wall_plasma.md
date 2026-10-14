# Resistive moving wall plasma

Adds prescribed scalar resistivity to the complete finite-mass, two-entropy
moving model. All ideal, Hall, current-stress/advection and species-pressure
terms remain. This is an opt-in bounded model; the collisionless solver and
its outstanding trajectory problem are unchanged.

For actual quadrature current J, add weak forcing `P.T W (eta J)` to the
original Ohm equation. The original inertia/Faraday solve therefore removes
magnetic plus current energy at `integral eta |J|^2 dV`. Let `Q=eta |J|^2`,
`S` be the actual fluid reconstruction, `V` its nodal volumes and `h_s` its
actual EOS internal-energy derivative. Deposit

    Cdot_s = f_s * (S.T W Q) / (V h_s),   f_e+f_i=1.

Then `sum V h_s Cdot_s = integral Q` by partition of unity. This is a local
Galerkin heat source paired with the resistive work, not a global energy
correction. Full grid power diagnostics include the heat; `ohmic_power`
reports the separate physical dissipation. Native versus continuously
reconstructed energy and finite-step drift still require measurement.

The high-order scalar adjoint is not positivity preserving: projected nodal
heat can be negative for poorly resolved Q even though quadrature Q is
nonnegative. No clipping, floor, filter or extra numerical diffusion is added.
The existing stage and endpoint positivity checks remain; smooth source
convergence and physical versus nodal heating must be evaluated.

`resistivity` is a finite nonnegative SI scalar. `electron_heat_fraction` is
required explicitly and is a prescribed model input, not a derived collision
partition. Temperature-dependent/tensor transport, electron-ion temperature
exchange, neutral collisions, ionization, conduction and radiation are absent.
Do not label this a complete Coulomb or laboratory transport closure.

`ResistiveTransport` plugs into the existing full RHS assembly so resistivity
participates in its single combined induction solve. Cached copies rebind the
actual density model. `ResistiveWallPlasma` keeps periodic axial and impermeable
radial boundaries, fixed poloidal wall flux and zero wall Btheta. It introduces
no vacuum, coil or formation model. Used by `examples.resistive_bessel_decay`.
