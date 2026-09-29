# Next moving-source coupling and spatial-accuracy checks

This is a derivation and implementation plan, not a completed coil-coupled
moving solver. Scalar resistive/moving results and their unresolved flow error
are frozen in RESISTIVE_MOVING_VERIFICATION.md. Hardware remains prescribed
ampere-turn waveforms; no circuit model is implied.

## External field at moving material coordinates

Let `a_ext(x,t)=psi_ext/r²` be the analytic circular-coil field plus any uniform
bias. The material total coefficient follows the existing flux-label rule;
the full-space induced problem prescribes

```
u_material = a_material(x, labels) - a_ext(x,t)
load = (2*pi/mu0) (K(x) u)_material
E_ind = (pi/mu0) u.T K(x) u
```

Stationary vacuum eliminates derivatives of the vacuum coefficients, but does
not eliminate the external field's dependence on moving coordinates. Thus

```
force = - partial_x E_ind at fixed coefficients
        - load dot partial_x a_material + load_i grad a_ext(x_i,t)
source_power = - sum_i load_i partial_t a_ext(x_i,t)
```

The first term includes the exterior mesh-motion transpose already implemented
by `MovingField`. The second is `MaterialFlux.pullback`. Omitting the third
would omit the nonuniform source force even if a uniform-bias budget passed.
Away from the axis `a_ext,r=(Bz_ext-2*a_ext)/r` and `a_ext,z=-Br_ext/r`.
On the regular axis these become zero and `B_axis'(z)/2`, respectively. Use
cancellation-free near-axis series rather than small radial division:

```
a_ext = B_axis/2 - r² B_axis''/16 + r⁴ B_axis''''/384 + ...
a_ext,r = -r B_axis''/8 + r³ B_axis''''/96 + ...
a_ext,z = B_axis'/2 - r² B_axis'''/16 + ...
```

At fixed geometry, the resistive update replaces the scalar bias increment by
the vector `Delta a_ext`: its forcing is `-M_full Delta a_ext`, and source work
is `-(2*pi/mu0) (K u_mid)_material dot Delta a_ext`. In the extended-time Strang
method, source time still advances only in the ideal component.

Evaluate physical B as induced finite-element B plus analytic external B.
Interpolating analytic coil fields into the induced stiffness would introduce
spurious external current. Initial material `a=a_ext` must give exactly zero
induced coefficients and recovered current. Coils must remain outside the
actual plasma; a singular filament crossing must fail atomically. A coil may
lie inside the magnetic matching region because its analytic field is separated
from the decaying induced field. Verify this explicitly.

## Independent exact force fixture

Inside a sphere of radius R, prescribe the induced coefficient
`a_self=A (1-3*s²/(5*R²))`, where `s²=r²+z²`. Outside use
`a_self=(2*A/5)(R/s)³`. Both a and its normal derivative are continuous, so
there is no surface current sheet. Ampere's law gives
`Jtheta=6*A*r/(mu0*R²)` inside and zero outside. Its entire exterior is a pure
dipole with moment `m=8*pi*A*R³/(5*mu0)`.

For a circular filament entirely outside this sphere, mutual-field reciprocity
gives the exact axial force `Fz=m B_axis'(0)`. This is exact for this spherical
current distribution, not a small-radius approximation. Independently integrate
`Fz=-integral_plasma Jtheta Br_external dV` in spherical coordinates. For coil
radius Rc, position Zc and ampere-turns Ieff,

```
B_axis(z) = mu0 Ieff Rc² / (2 (Rc²+(z-Zc)²)^(3/2))
B_axis'(z) = -3 mu0 Ieff Rc² (z-Zc) / (2 (Rc²+(z-Zc)²)^(5/2))
```

A direct 48×48 Gauss volume check using the independently existing loop-field
formula at R=1 m and A=sqrt(mu0) T gives:

| Rc [m] | Zc [m] | Ieff [A turn] | Exact Fz [N] | Relative integral difference |
|---:|---:|---:|---:|---:|
| 1.6 | 0.8 | 1200 | 1.133971916017153 | 2.7e-15 |
| 2.0 | -0.5 | -700 | 0.3177764850317322 | 2.2e-15 |
| 0.4 | 1.7 | 950 | 0.1344339493450913 | 3.6e-15 |

The 24/48/96-order integral check converges to these values. Turn this derived
reference into a durable independent benchmark and tests before implementation.
Initialize total `a=a_self+a_ext`, then compare net axial virtual force under
mesh, matching and quadrature refinement. Also perturb the complete potential
under arbitrary radial/axial motions to check every geometric derivative. This
fixture is instantaneous force verification, not a freely equilibrated plasma.

## Resolve the existing spatial error separately

The current D=0.1 driven nr=32 speed maximum lies one material shell inside the
equatorial surface, while nr=8 peaks on the surface. D=1 and unforced resistive
relaxation peak on the equatorial surface. Localizing error there does not
identify its cause. Inspect radial velocity/pressure/current profiles and
discrete force contributions, surface heat, axis effects and material modes.
Refine time at the finest mesh if temporal instability remains a concern.

Initial FRC pressure vanishes at the surface but current does not. Ohmic
heating changes that surface thermodynamics immediately. A resolved heated
free-surface expansion is a distinct requirement from fixed-sphere diffusion;
the current measurements do not validate its continuum limit. Develop an
independent gas/vacuum heating or expansion reference and check the boundary
traction/entropy treatment. Do not suppress the error with a pressure floor,
undisclosed artificial resistance or a change of physical D. Keep formation
gated until current, flow, boundary and transport assumptions are resolved.
