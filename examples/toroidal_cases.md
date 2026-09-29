# Nonlinear toroidal measurement cases

The 33-case protocol contains nine independent rotation/pinch equilibria,
eight contact/wave axial refinements, three nonlinear contact controls, four
wave amplitude/mode/order controls, five conducting-Bessel diffusion refinements
through nr=128, and four separate time caps. Complete physical/model/mesh inputs
are explicit. The reference boundary is never substituted for free vacuum.

`build` creates a fresh solver. `reference` independently constructs the expected
conservative averages or asymptotic field at the requested physical time.
Contrasts are finite amplitude for the exact contact; weak waves and diffusion
carry their explicitly stated nonlinear feedback error. All use SI mathematical
fixtures, not newly chosen laboratory parameters. No result is FRC formation.
The requested gamma is explicitly passed through every fixture energy
constructor, including torsional and diffusion references.
