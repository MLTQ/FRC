# Analytic spherical FRC in genuine magnetic vacuum

This independent continuum fixture supplies a target for future free-boundary
verification. It is **prescribed static equilibrium**, not evolution from coils,
not a density/transport closure and not evidence of stability or confinement.
The earlier dynamic solvers retain their declared fixed chamber boundaries.
The new ideal material reference evolves this fixture with an actual moving
surface; its drift/refinement is reported in MOVING_PLASMA_VERIFICATION.md.

The poloidal Hill-vortex form is used as an FRC reference in
[Glasser and Cohen (2022), section III.A](https://w3.pppl.gov/ppst/docs/glasser2022rsi.pdf).
The spherical matching construction is discussed by
[Keller and Cheviakov (2022)](https://arxiv.org/abs/2204.02192).
Below we independently normalize the spherical solution to this repository's
SI convention and an imposed positive or negative axial bias B0 at infinity.
No stability conclusion is imported from these references.

## Fields and pressure

Let s²=r²+z², separatrix radius a>0 and C=-3 B0/4. With
Br=-psi_z/r and Bz=psi_r/r, prescribe

```
psi_inside  = C r² (1-s²/a²),                         s <= a
psi_outside = (B0/2) r² (1-a³/s³),                    s > a
p_inside    = (10 C/(mu0 a²)) psi_inside >= 0
p_outside   = 0
Jtheta_inside = 10 C r/(mu0 a²); Jtheta_outside = 0.
```

Inside, Delta*psi=-10 C r²/a². Thus Jtheta=-Delta*psi/(mu0 r)
and J×B=(Jtheta/r) grad psi=grad p. The exterior is uniform axial
bias plus a decaying dipole, with zero curl B. There is no toroidal field,
flow or material outside the sphere. The regular axis limits are evaluated
directly, without dividing by r.

Both pieces have psi=0 at s=a. Their normal derivatives agree because
`-2 C/a = 3 B0/(2a)`. Therefore B is continuous and there is no surface
current sheet. On the sphere, in polar angle theta,

```
Br = -(3/2) B0 sin(theta) cos(theta)
Bz =  (3/2) B0 sin(theta)²
B_normal = 0; p = 0.
```

The full magnetic traction is continuous; the normal total pressure agrees
without exterior gas pressure. Volume current has a finite jump, so derivative
tests must distinguish smooth interior/exterior points from the interface.
Zero physical resistivity is required for this unsustained static equilibrium.
Finite eta would produce induction and heat; it is not silently balanced by an
invented electric drive.

## Independent topology and integral targets

The center-axis field is Bz=-3 B0/2. The magnetic null/O-ring lies at
`r=a/sqrt(2), z=0`, with `psi_O=-3 B0 a²/16`. The spherical separatrix
ends at axis nulls z=±a. In our flux-per-radian convention, the magnitude
of closed poloidal flux is `2pi |psi_O|=3pi |B0| a²/8`.

Writing S=pi B0² a³/mu0, direct polynomial and dipole integration gives

| Quantity | Exact continuum value |
|---|---:|
| Separatrix volume | 4 pi a³ / 3 |
| Integral of pressure | 6 S / 7 |
| Gas internal energy for gamma>1 | 6 S / [7 (gamma-1)] |
| Total magnetic energy inside the sphere | 3 S / 7 |
| Induced magnetic energy over all space | 10 S / 7 |
| Induced magnetic energy outside the sphere | S / 3 |
| Toroidal current integrated over a poloidal half-plane | -5 B0 a / mu0 |

The uniform-bias self energy over infinite space diverges and is excluded from
the induced-energy target. These are field/plasma integrals, not coil-circuit
energies. Independent Gauss quadrature uses x=a/s to cover the entire vacuum
without imposing a finite outer boundary or fitting its missing energy.

## How the reference is used

The independent formulas in `examples/spherical_frc_reference.py` are separate
from numerical evolution. Tests verify differential force balance, vacuum
current, axis/flux signs, interface continuity and full-space energy. The ideal
material prototype now initializes these fields through Q1 magnetic values and
quadrature-integrated cell pressure, and measures force, shape, flux and energy
drift under spatial/time refinement. Its deliberately uniform verification
density does not select a physical species or transport regime. Perturbation dynamics
require their own reference or convergence tests; the static target is not
a proof of stability.

Density, species and electron/ion temperatures are deliberately unspecified:
magnetostatic force balance alone does not determine them. Pressure vanishes
on the sphere and symmetry axis. A physical thermodynamic/kinetic closure near
those regions must be assessed instead of filling vacuum with an atmosphere or
adding a temperature/density floor. No laboratory operating regime follows from
the unit-scale verification constants.

## Measured reference checks

`python -m examples.verify_spherical_frc` records independent continuum energy
quadrature and applies the existing topology detector to this exact field.
At 8/16/24 Gauss points, all four energy/pressure integrals agree with the
closed-form targets within 2.11e-15 relative error. This is integration of a
known formula, not conservation by a time-dependent solver.

| Radial intervals | Separatrix-volume relative error | O-rings found | Axis nulls found |
|---:|---:|---:|---:|
| 24 | 0.00219889 | 1 | 2 |
| 48 | 0.000606522 | 1 | 2 |
| 96 | 0.000157794 | 1 | 2 |
| 192 | 0.0000374111 | 1 | 2 |

The detector's axis-attached separatrix estimate converges approximately at
second order. Its largest sampled closed contour has smaller volume and must
remain a lower bound, not be relabeled as the separatrix. Noise, degenerate
critical points and subgrid features remain separate verification requirements.
Complete reference parameters, detector outputs and source/environment hashes
are saved in `outputs/free_boundary_reference/metrics.json`.

![Prescribed spherical reference](../outputs/free_boundary_reference/spherical_frc.png)
