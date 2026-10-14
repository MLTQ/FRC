# Free plasma expanding into its magnetic vacuum

Complete simultaneous evolution of the full material fluid/current/species state
and penetrating exterior flux through the moving side and both ends. Exterior
geometry and magnetic minimization are recomputed at every RK stage. The exterior
has exactly zero matter/conductivity. This is an unforced initialized field, not
coil formation or a validated laboratory plasma.

Initial psi=r²*f(s²), s²=r²+z². For a=.7R, u=s²/a²,
f=C/a³*(35-42u+15u²)/8 inside s<=a and f=C/s³ outside. Value and first two
derivatives match at s=a. Thus Jtheta=52.5*C*r*(1-u)/(mu0*a⁵) inside and zero
outside its spherical support, with an analytic dipole vacuum. C=Bcentral*a³/8.75.
There is no singular current sheet in these continuum initial data. The entire
plasma boundary lies outside the current support. Spatial projection errors
and weak boundary-current behavior remain measured limitations.

The independent initial acceleration includes analytic pressure, J cross B and
the radial curvature of the pure toroidal current stress. Vacuum initial energy
uses an independent Green-identity integral of the exact dipole. The continuum
all-space virial I''=4E-2Wb includes exterior magnetic energy. Its observer does
not repair evolution; finite exterior approximation can leave virial/translation
errors even when paired discrete power cancels exactly.

Every accepted state and source closure is retained. Prescribed resistance,
upwind entropy diffusion and conditional Coulomb exchange stay explicit. The
inherited upwind current time sensitivity remains unresolved. All physical
mass/heat/current terms and residual limits remain; no energy/state repair,
atmosphere, forced surface motion or prescribed shape evolution is added.

The first completed run archive used the temporary module name
`material_vacuum_expansion`; its frozen source is preserved. This distinct entry
point preserves the earlier nonmagnetic affine-expansion driver at that name.
