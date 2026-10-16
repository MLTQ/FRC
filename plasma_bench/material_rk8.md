# Fixed higher-order central plasma advance

`dop853_step(model, state, dt, time)` applies SciPy's twelve-stage DOP853
eighth-order tableau to the complete existing material derivative. Every stage
evolves coordinates, momentum, rotation, both entropies and induced magnetic
coefficients with actual prescribed-source time, collision heat and forces.
State addition retains the existing axis/boundary constraints; invalid physical
stages reject through the same model. No field, particle, entropy or energy
projection, floor, repair or extra exterior material is introduced.

The explicit equations are k_i=F(t+c_i*dt,y+dt*sum_j A_ij*k_j), followed by
y_new=y+dt*sum_i B_i*k_i. A/B/C are taken directly from the installed exported
`scipy.integrate.DOP853` class (SciPy 1.18.1 for the initial runs). This uses
only the principal tableau: no adaptive tolerance, embedded estimate or dense
interpolation is claimed. The timestep remains an explicit input. Stability
and resolution of fast current modes still require actual trajectory evidence.

[SciPy DOP853 documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.DOP853.html)
describes the principal method and links Hairer's original implementation.
Library version is retained by the existing source provenance recorder.

The caller is `EnergyDiffusionMaterialCoilCurrent.central_step` when
`central_integrator="dop853_fixed"`; omitted/`"rk4"` retains the prior method.
Conservative numerical thermal mixing still surrounds the central step, so
the full split evolution is at most second order even for smooth solutions.
Changing the central tableau does not establish full spatial accuracy, physical
thermal transport, vacuum-edge validity or FRC formation.
