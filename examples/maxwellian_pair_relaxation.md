# Complete homogeneous finite-drift collision relaxation

Run `python -m examples.maxwellian_pair_relaxation OUTPUT`. Two declared
collision-only problems evolve both three-component velocities and both scalar
temperatures through 100 initial small-drift relaxation times, with 100/200/400
RK4 steps. Equal masses expose comparable species heating; physical electron/argon
masses expose drag and equilibration timescale separation. Both start at relative
speed three times the combined thermal speed, with a common nonzero boost.
Density is fixed and the Coulomb logarithm is explicitly 10.

The independent reference integrates the binary momentum/energy moments directly
in relative-velocity spherical coordinates with 96 radial and 64 half-angle Gauss
points. The other three velocity integrals use the conditional Gaussian pair
centre velocity. The tail ends 12 thermal speeds beyond the relative drift;
opposite angles are paired using expm1 to preserve tiny drift. No production
collision function, erf or incomplete gamma is used by this reference. DOP853
advances these velocity-integral rates with nondimensional state/time tolerances
rtol=2e-11, atol=2e-13. Comparison uses every accepted RK4 state and actual physical
energy, momentum, entropy, drift and temperatures; no energy correction is made.

Total initial momentum fixes common eventual velocity. Subtracting its kinetic
energy from initial total energy fixes the eventual common temperature. That is
an independent equilibrium prediction, not an imposed state. Physical Ar masses
can leave temperatures far from equilibration over this drag-timescale window.

This is a complete collisional-source trajectory, not the full electromagnetic
moving-plasma solver or an experimental formation baseline. Fields and spatial
transport are deliberately absent from the stated homogeneous operator problem.
Maxwellian shapes are assumed throughout; distribution relaxation, anisotropy,
screening/cutoff choice, kinetic instability and collision magnetization are not
established by reproducing instantaneous moments. Existing coil runs are untouched.

The output freezes source/companions, saves all accepted histories and independently
integrated references, and emits one scientific PNG/PDF and concise measurements.
No unit tests, pass/fail suite or recursive verification framework are introduced.
