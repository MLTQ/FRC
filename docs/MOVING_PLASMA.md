# Experimental ideal plasma with a moving material surface

This reference evolves an actual two-dimensional axisymmetric plasma boundary
against a massless magnetic vacuum extending to infinity. It is an **ideal,
adiabatic, smooth-flow prototype**, separate from the fixed-chamber resistive
solvers. It has no exterior gas, pressure/density floor, numerical viscosity,
numerical resistivity or energy repair. The Stage 2 formation entry remains gated.

The continuum starting point is ideal material MHD: fixed parcel mass, advected
entropy, advected poloidal flux and magnetic virtual work. Material variational
MHD is discussed by [Zhou et al. (2014)](https://arxiv.org/abs/1408.1346).
The present Q1/P0 axisymmetric discretization is **not their discrete exterior
calculus scheme** and does not inherit their topology/momentum guarantees.
The usual moving plasma–vacuum problem has material interface motion, continuous
total traction and a curl-free exterior field; see
[Secchi and Trakhinin (2013)](https://arxiv.org/abs/1301.5238).
Their non-collinearity stability hypothesis is not established for this
axisymmetric FRC fixture. No well-posedness or stability theorem is imported.

## Governing quantities and regular magnetic variable

Use SI and flux per radian, so the physical poloidal flux is 2 pi psi:

```
B_r = -psi_z/r,   B_z = psi_r/r,   B_theta = 0.
D rho/Dt = -rho div v,
rho Dv/Dt = -grad p + J cross B,
D(p/rho^gamma)/Dt = 0,   D psi/Dt = 0,
J = curl B / mu0,   curl B = 0 in vacuum.
```

Set `a=psi/r²`, with a regular finite axis limit. Then
`Br=-r a_z`, `Bz=2a+r a_r`, and vacuum satisfies
`div(r³ grad a)=0`. This removes the singular `1/r` coefficient of a direct
psi discretization. Continuity of `r²a` across element faces makes the normal
magnetic flux continuous in the weak sense; cellwise div B is zero analytically.

For any regular decaying induced field, expansion of B² gives

```
r [(2a+r a_r)² + (r a_z)²]
    = r³ |grad a|² + partial_r(2r²a²).
E_ind = (pi/mu0) integral_halfplane r³ |grad a_ind|² dr dz.
```

The divergence term vanishes at the axis and infinity, and cancels across
internal interfaces. Thus the last formula equals the physical full-space
induced B²/(2mu0) energy. The implementation independently integrates both
forms as a quadrature check; the agreement is not itself a spatial-accuracy test.
Uniform external-field self energy is infinite and is excluded. The finite
response energy is material kinetic plus internal plus induced magnetic energy.
It is not a coil-circuit or delivered hardware-energy estimate.

## Material mesh, thermodynamics and advection

`MovingDomain` begins with concentric radial shells and angular sectors in the
meridional half-plane. Q1 bilinear cells revolve about r=0; the shared center
has an intentionally collapsed inner edge. The material surface is a polygon
approximating a sphere. Every material vertex moves in r and z, except that
axis nodes remain at r=0. There are no exterior fluid variables.

Cell masses `m_c=integral rho_initial dV` and positive lumped nodal inertial
masses `M_i=integral rho_initial N_i dV` are fixed. Pressure is cellwise constant:

```
kappa_c = p_initial_cell_average V_initial^gamma,
rho_c = m_c/V_c,  p_c = kappa_c/V_c^gamma,
U_gas = sum_c kappa_c V_c^(1-gamma)/(gamma-1).
```

No mass leaves a material cell. This P0 entropy approximation needs refinement
where the prescribed initial pressure varies. It is not a two-temperature,
ionization or collisional transport closure.

At non-axis material vertices, the exact nodal flux label gives
`a_i=a_initial_i (r_initial_i/r_i)²`. At the axis psi=0 alone does not specify
Bz. The axis a value is divided by the square of an estimated neighboring
radial stretch: the nearest off-axis node at each pole, and a sin²-weighted
average on the first shell at the center. This is exact for a uniform radial
affine deformation, but is an explicit spatial approximation otherwise.
The force includes its exact derivative. Nodal flux invariance does **not**
prove frozen flux on every interpolated moving surface or preservation of all
topology. Resistive topology change and remapping are absent.

## Finite magnetic mesh and exact infinite exterior

Vacuum nodes between the moving material boundary and a fixed sphere of radius
R follow the linear extension `(R-s)/(R-a0)` of the boundary displacement at
the same initial angle. This is a computational mesh rule, not vacuum inertia.
The last magnetic annulus maps exactly to the sphere using the analytic arc;
its field basis blends inner Q1 traces with exact angular modes on the sphere.
Continuity holds over entire faces, not merely at vertices.

With mu=z/s and `norm_l=2l(l+1)/(2l+1)`, l>=1, use

```
Phi_l(mu) = P_l'(mu)/sqrt(norm_l),
a_l(s,mu) = c_l (R/s)^(l+2) Phi_l(mu).
integral_-1^1 (1-mu²) Phi_l Phi_k dmu = delta_lk.
```

Separation of the weighted vacuum equation gives decay exponent l+2. Integrating
the radial and angular energies gives the exterior weighted diagonal
`(l+2) R³`; direct physical B² energy outside R gives `l R³`. Their difference
is the boundary-divergence term above and cancels with the finite-domain term.
Both receive the prefactor pi/mu0. Conflating these two diagonals is incorrect.

The assembled symmetric matrix K contains finite-element stiffness and the
analytic exterior diagonal. Prescribe `a_ind=a_material-B_external(t)/2` at
material nodes; solve `K_vv a_v=-K_vp a_p` for all vacuum nodes and modes.
There is no exterior electrical conductivity. The vacuum residual measures
algebraic stationarity; matching radius, angular/radial spacing, number of
modes and quadrature must be independently varied to assess accuracy.

## Magnetic and pressure forces

The potential is `U=U_gas+(pi/mu0) a_ind^T K(x) a_ind`. Its negative derivative
with respect to material positions supplies the force. Vacuum stationarity
allows the envelope theorem: the derivative of the solved vacuum coefficients
does not contribute, but all vacuum geometry derivatives do.

For a geometric shape function N and element magnetic gradient g, the analytic
fixed-coefficient derivative of the weighted magnetic integral is

```
integral [r³ (|g|² grad N - 2(g dot grad N)g)
          + 3r² N |g|² e_r] dr dz.
```

The transpose of the declared exterior mesh-motion map returns vacuum geometric
terms to material vertices. The additional material coefficient derivative is
`load_i da_i/dx`, where `load=(2pi/mu0) K a_ind`; off axis `da_i/dr_i=-2a_i/r_i`.
The regular-axis stretch rule supplies its own chain derivative. Gas force is
`sum_c p_c grad_x V_c`. These analytic derivatives are checked against an
independent symmetric perturbation of the complete potential and against a
closed-form sphere energy derivative. They are not adjusted to hold a fixture
in equilibrium. The radial axis reaction is constrained and performs no work.

At a tangent-field material interface, continuum variation recovers the
magnetic and gas traction balance. Finite mesh resolution does not enforce
that balance exactly. A uniform magnetized sphere with a tangential field jump
and nonzero normal B is an **electromagnetic virtual-force check only**, not a
valid free-gas equilibrium. The separately derived spherical FRC has p=0 and
Bn=0 at its surface with continuous B, so it supplies the interface target.

The continuum energy in uniform bias is invariant under axial translation.
The finite matching mesh generally is not, producing a measured spurious net
axial force. Exact magnetic axial-momentum conservation is therefore not claimed.

## Time stepping and prescribed-source work

For fixed masses use velocity Verlet:

```
v_half = v_n + (dt/2) F(x_n,t_n)/M,
x_new = x_n + dt v_half,
v_new = v_half + (dt/2) F(x_new,t_new)/M.
```

At fixed material coordinates the source contribution is
`P_source=-(Bdot_external/2) sum_material load_i`. Integrate it by the
trapezoid rule. The response diagnostic is
`K+U_gas+E_ind - E_initial - W_source`; finite-step error is measured, not repaired.
Second-order time behavior is tested separately from spatial equilibrium error.

The step cap uses cell meridional width `2 area/perimeter`, sound plus magnetic
fast-speed bound and relative node velocities, with CFL<=0.5. Arbitrary dynamic
sources require an explicit positive additional timestep cap. `step(dt)` is
also available for controlled time studies; callers then own the chosen dt.
Geometry, field and velocities must validate before state, clock or cumulative
work are committed. Corner and Gauss Jacobians are checked, allowing only the
intentional center collapse. This is not a global mesh-injectivity theorem.

## Exact expanding-sphere reference

For B=0, gamma=5/3, initial rho=1 and p=1-s² inside a unit sphere, take material
position `x=lambda(t) X`. Then rho=lambda^-3 and p=(1-|X|²) lambda^-5. The
momentum equation gives `lambda''=2 lambda^-3`; with lambda(0)=1, lambda'(0)=0:

```
lambda=sqrt(1+2t²),   v=(2t/lambda) X.
```

The surface pressure is zero and no external fluid is required. This verifies
actual material motion against a continuum solution, not just a budget identity.
The spherical FRC reference instead has zero exact acceleration and motion;
its measured drift is a numerical error. Both use deliberately unit-scale
mathematical density/pressure choices, not laboratory argon parameters.

## Current boundaries of validity

The prototype lacks resistive induction/heating, spatially varying coil fields,
shocks and entropy production, mesh remapping, arbitrary chamber geometry,
neutrals/ionization/radiation and Hall/two-fluid physics. Q1/P0 material meshes
can have low-energy distortion modes and will eventually distort under strong
flow. There is no hidden damping that cures them. No large-deformation,
long-time stability, physical FRC formation or confinement claim follows from
the current smooth short-time references. A nonideal extension must retain
vacuum/source work and add independently verified heat/entropy accounting.
