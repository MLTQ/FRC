# Exact polynomial subdivision of diagnostic directions

`transfer(source, target, value)` evaluates the original conservative direction
on nested target fluid nodes with local root-product Lagrange polynomials in
x=r²/2. Child ownership determines both parent traces. Scalar components use
degree five; transverse regular momenta use the dependent-axis degree-four
polynomial in the first parent element. New inactive-axis slots are omitted.
Parent endpoint coefficients are copied exactly as prescribed trace data.

No new direction normalization, least-squares fit, physical state projection,
filter or mass change. The source/target axial nodes and degree must match.
The source is an archived diagnostic direction, never an evolved state.
`parent_faces` selects the two owning-parent endpoint arrays for verification.

At the radially constant base this preserves differential physical velocities.
Finite states base±h*direction can have different interpolated nodal q/n fields
because that representation is nonlinear. The separate end-to-end transfer audit
must measure those differences for every amplitude/sign before interpreting
coupled rates; exact conservative subdivision is not claimed to preserve every
finite-amplitude primitive interpolation identically.
