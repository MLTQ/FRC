# Magnetic vacuum to infinity

## Why this boundary is needed

The fixed zero-induced-flux box changed the coil-pulse source work by about 19%
when its extent doubled. It accurately solved the wrong finite exterior for a
free-space experiment. `vacuum_boundary="multipole"` now joins the finite field
grid to an analytic decaying vacuum expansion. Its magnetic energy beyond the
grid is retained in the dynamic work identity. The old fixed option remains a
comparison, and the fixed **fluid chamber** is unchanged by this improvement.

## Exterior representation

Let R=sqrt(r²+z²), x=z/R, a=min(box radius, box half-length). For l=1,...,L use

```
F_l(r,z) = (a/R)^l (1-x²) P_l'(x) / sqrt(l(l+1)),
psi_exterior = sum_l c_l F_l.
```

P_l is a Legendre polynomial. Each mode solves Delta*psi=0, vanishes as r² on
the axis and decays at infinity. Dipole l=1 gives a r²/(sqrt(2) R³). Odd and even
modes are retained, so axial asymmetry is represented. Analytic derivatives give
Br=-psi_z/r and Bz=psi_r/r, with explicit regular limits for exterior axis probes.
The expansion is only for **induced** fields; prescribed analytic coils/bias are
evaluated independently and need not decay or be inside the matching surface.

All induced current lies inside the fluid chamber. The matching rectangle must
enclose the sphere R<=sqrt(chamber radius²+chamber half-length²); the constructor
rejects geometries that do not. This ensures that the exterior expansion is
outside every possible induced source. It can require a relatively wide field
grid for a long chamber. Multiple expansion centers or other exterior methods
may be more efficient in that geometry, but are not silently substituted here.

## Variational coupling and energy

The interior flux stencil has a symmetric Dirichlet energy. Partition all
non-axis nodal values into interior I and outer boundary b. Boundary radial/axial
tangential edges use half-cell weights; this gives the full interior energy

```
u_I^T K_II u_I + 2 u_I^T K_Ib u_b + u_b^T K_bb u_b.
```

Set boundary nodal values u_b=F_b c. The exact continuum exterior energy of the
retained modes is c^T T c. Green's identity gives

```
T_lm = integral_exterior grad(F_l) dot grad(F_m) / r dr dz
     = - integral_boundary F_l (dF_m/dn) / r ds,
```

where n points **out of the finite box**, opposite to the exterior domain's inner
normal. The axis and infinity contribute zero. Integration over the radial side
and both end faces uses Gauss quadrature; reciprocity and positive eigenvalues
are checked. An independent dipole volume integral verifies the sign and scale.

With q=(u_I,c), the full magnetic energy is W=pi/mu0 q^T K_open q, where

```
K_open = [ K_II          K_Ib F_b                 ]
         [ F_b^T K_bI    F_b^T K_bb F_b + T       ].
```

This matrix is symmetric. Exterior nodal rows and modal rows impose stationary
energy (K_open q)_vac=0. Modal rows have no fluid mass, resistivity, current
density or velocity interpolation. The existing midpoint/advection adjoint pair
therefore preserves delta(W+E_gas)=applied-field work with **the energy outside
the grid included**. Modal stationarity residuals are reported separately from
physical vacuum current density. No sheet current or artificial outer wall is
prescribed on the matching surface.

## Approximation and verification contract

The exterior satisfies the vacuum differential equation exactly within the
retained basis. A finite L truncates that basis. Interior edges interpolate its
nodal trace, whereas the exterior trace is analytic; this introduces spatial
consistency error at the interface. Weak matching in retained modes does not
promise pointwise derivative continuity on a coarse mesh. Independent normal-
derivative jumps, mesh refinement, mode refinement and matching-surface movement
must be measured. Energy conservation alone cannot establish any of these.

`open_benchmarks` excites both parities using a compact displaced toroidal current.
Its reference integrates exact circular-loop fields independently of K_open.
The checks include the field inside the grid, far outside it and on the exterior
axis. Source quadrature error, spatial error and mode truncation are recorded
separately. Dynamic pulse comparisons retain identical gas, source and time-step
parameters while changing only the exterior boundary choice/surface/mode count.

Related open-domain magnetostatic work includes
[Bonotto, Abate and Bettini's axisymmetric FEM–BEM coupling](https://www.research.unipd.it/handle/11577/3478668)
and [Birk, Reichel and Schröder's exterior energy reduction](https://doi.org/10.1016/j.cma.2022.115362).
Those works motivate coupling a finite field domain to its unbounded exterior.
This nodal/multipole scheme is independently derived; it is not an implementation
of their finite elements or a claim to inherit their convergence results.

An open **magnetic** boundary does not yet give a moving plasma surface, gas
outflow, realistic quartz/end physics, physical transport, magnetic-shock fidelity
or FRC stability/confinement. Those remain separate requirements.
