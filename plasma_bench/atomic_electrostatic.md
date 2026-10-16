# Spatial charge, longitudinal field and compatible work

`AtomicElectrostatic` extends the actual species Euler transport with a 1D
longitudinal electrostatic field. Each cell retains every heavy/electron number,
three momenta and total translational energy. Reactions/collisions may use the
inherited source splitting, but the initial wave study disables them explicitly.
This is a closed reflecting verification interval, not a quartz-wall, sheath,
physical open-end or FRC free-boundary model. The formation exterior is unchanged.

The electric field is a Gauss-law unknown determined from the actual charge:

    E[j+1]-E[j] = dx rho[j]/epsilon0,
    rho[j] = e Nref sum_s Z_s f_s[j].

The left normal field is a fixed explicit input (zero in the wave cases); the
right field follows the actual total charge. No mean charge is subtracted and
no neutralizing background is inserted. This constraint solve defines the
electrostatic model; it does not project an independently evolved Maxwell field
or alter any species state. Reflecting lab-stationary walls impose zero particle
and lab-energy flux; their normal pressure stress remains. Fixed-frame energy
flux includes the transformed wall pressure work.

Let F_s[j] be the ACTUAL numerical particle flux, including the Euler numerical
diffusion. The face charge current is J[j]=e Nref sum_s Z_s F_s[j]. Differentiating
the same Gauss equation and continuity gives dE[j]/dt=-J[j]/epsilon0, since the
left wall has J=0. Thus field evolution and particle transport use one current,
not separately reconstructed inconsistent moment and density currents.

Cell momentum source is q_s n_s (E_left+E_right)/2. Species lab-energy source in
eV per reference nucleus is

    W_s = Z_s (E_left F_s,left + E_right F_s,right)/2.

Its sum over cells cancels the derivative of the face-lumped field energy
epsilon0/2 * trapezoid(E²,dx). This is a semidiscrete source derived from the
numerical continuity law, not an after-the-fact global energy correction.
Species fixed-frame work subtracts U dot dP. Field energy is computed separately;
SSP-RK2 time integration need not preserve the quadratic total exactly.

This work contains numerical flux contributions; it must not be called solely
physical Joule heating. It consistently approximates q n v E as the mesh is
refined. The second-order face-lumped energy differs from the exact integral of
the piecewise-linear Gauss field. `totals` reports both budgets, plus the actual
Gauss residual. `last_diagnostic` measures the semidiscrete Ampere and work
residuals independently from the evolved density rates. Neither changes a state.

The timestep resolves physical electron acoustic and local plasma-frequency
scales via explicit CFL and an explicitly supplied maximum omega_p*dt. No
artificial electron mass, density floor, clipped charge or field-energy repair.
HLLC/HLL/minmod diffusion is inherited and disclosed. Interior vacuum tails from
that Eulerian method remain unsuitable for replacing the tracked formation edge.

The general finite-volume issue that field work must follow numerical matter
fluxes is also discussed for self-gravity by
[Springel (2010)](https://doi.org/10.1111/j.1365-2966.2009.15715.x).
The electric equations and boundary identities above are the explicit derivation
used here; gravitational sign conventions are not imported. Transverse Maxwell
curls, magnetic forces, radial geometry, open particle loss and a compatible
moving interface remain necessary for the full sandbox.

`AtomicMaxwell` now extends this longitudinal operator with transverse fields;
see atomic_maxwell.md. It depends on the same Gauss/current/work and reflecting
flux contracts. Completed longitudinal measurements are in
`docs/ATOMIC_SPATIAL_FIELDS.md`; neither component is a formation-edge replacement.
