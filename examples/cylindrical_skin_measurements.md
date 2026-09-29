# Kinematic cylindrical skin verification

`RadialSkin` evolves the invariant-z magnetic flux of a stationary, uniform
Ohmic cylinder. It uses the actual production `vacuum_flux.flux_components`
and shared-entity fallback for the `standard` and `strip` cases. The latter
forces the last radial cell to first order. Eight identical axial cells provide
the production stencil; only the central plane is retained after each rate
evaluation, enforcing exact axial invariance. End and corner dynamics are not
part of this reduced test. No finite-cylinder or full MHD claim is made.

Gas is held at rho=1, p=1, v=0 to declare the production LLF speed. The tiny
field amplitude makes the speed nearly constant, but the benchmark is a
stationary conductor by definition, not an evolved approximate plasma. The
exterior field is prescribed uniform Bz(t), with no exterior current. The
radial stiffness current uses centered B differences and a half-cell boundary
dual width; the surface flux evolves, and Bz(R,t) supplies natural boundary
data. The exterior induced B is exactly zero for this infinite-cylinder model.

`physical` omits all UCT and interface numerical diffusion and independently
checks the radial physical diffusion baseline. Other methods retain the
production numerical wall eta=a*h*mu0/4 and LLF reconstruction. SSP RK3 evolves
the induction equation with a conservative explicit diffusion/advection cap.
No hydrodynamic or total-energy update is performed in this reduced experiment.

`measure` initializes the independent exact periodic Bessel solution and
compares eight times over one period. Exact B is averaged by analytic flux
differences over the cylindrical cell. Current errors use the continuum
derivative, including the one-sided interior surface current. All errors use
fixed applied-field/skin scales, never a near-zero instantaneous field norm.
`run` measures 16/32/64/128 and half timestep at 64, saves complete profiles,
inputs, method hashes and engine provenance. This tests spatial skin diffusion
and numerical dissipation; the paired full coil-pulse study tests gas response.
