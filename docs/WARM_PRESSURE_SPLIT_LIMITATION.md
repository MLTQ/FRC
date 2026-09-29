# Pressure-only splitting: a retained failure and a continuum warning

The new conservative warm integration is not yet a generally validated solver.
Fine pressure-only evolution rejects a negative ion-energy partition even while
electron pressure is positive, the nonlinear equations converge and global
energy is conserved. Smaller time steps retain the failure. An independent
local continuum symbol also has nonreal characteristic speeds for the same
nonbarotropic reference. These are distinct pieces of evidence; the symbol
counterexample does not by itself prove the cause of the observed axis error.

## Actual rejected evolution

The original ElectronReference has R=.035 m, L=.24 m, nref=1e20 m^-3,
Te scale=5 eV, density contrast .4, temperature contrast .3, toroidal slope
.05 T/m and gamma=5/3. Density/ions are fixed in this subproblem. The initial
uniform ion pressure is 16 Pa. The pressure midpoint evolves U and g; total
energy advances only by the conservative pressure Poynting/enthalpy flux.
The requested interval is 1 microsecond. No floor or energy repair is applied.

At nr64, doubling planned steps 80→160→320→640 changes the maximum local
ion-pressure error from .0920654 to .0915853 to .0914659 to .0914361 Pa.
The error remains as time error shrinks. Volume RMS remains about .0011633 Pa,
so this local problem would be easy to miss in an integrated norm.

At nr128/nz256 the following trials are rejected:

| Planned steps | Accepted steps | Last accepted time [µs] | Rejected time [µs] | Trial minimum pi [Pa] |
|---:|---:|---:|---:|---:|
| 320 | 303 | .946875 | .950000 | -.495657 |
| 640 | 607 | .9484375 | .950000 | -.612730 |
| 1280 | 1214 | .9484375 | .94921875 | -.315896 |

All three minima occur in radial/axial cell [0,100]. Trial electron pressure
is above 51 Pa and nonlinear residuals are 4.42e-11, 2.26e-11 and 1.05e-11,
below the declared 1e-10 tolerance. Maximum accepted-history relative total
energy drift is 2.08e-16. Global ion energy also remains constant: local
numerical redistribution creates the invalid cell, not net collisional heat.

The original supplemental driver raises before archiving this failed case.
`warm_pressure_failure` reproduces it with prospective partition validation,
preserving initial, last accepted and rejected arrays under separate keys.
The trajectory `completed` flag is false; `experiment_complete` only means
the failure observation finished. All histories, source/configuration/archive
hashes and independent rejected-pressure reconstruction pass their audit.

![Rejected pressure-only trials](../outputs/warm_pressure_failure/pressure_split_rejections.png)

## Independent frozen principal symbol

Take a smooth local point away from the axis, freeze positive n, pe and their
gradients, and retain only derivatives of perturbations. Let a=grad log n,
c=grad pe-gamma pe a, and let (kr,kz) be a **unit** wavevector direction.
The isolated continuum pressure equations are

    Btheta_dot = (nr pe_z - nz pe_r)/(e n²)
    pe_dot = J/(n e) dot (grad pe - gamma pe grad log n).

Here nr,nz denote partial derivatives, not grid counts. Linearizing gives

    delta(Btheta,pe)_dot = i |k| M delta(Btheta,pe)
    M = [[0,A],[C,D]]
    A = (ar kz - az kr)/(n e)
    C = (kr cz - kz cr)/(mu0 n e)
    D = (Jr0 kr + Jz0 kz)/(n e).

The background-current advection D is included. Cylindrical curvature and
terms without perturbation derivatives are lower order in this local limit.
M has real characteristic speeds only if D²+4AC is nonnegative. The units of
Btheta and pe differ, but this discriminant is invariant under their rescaling
and has units m²/s². If it is negative, one eigenvalue of i|k|M has positive
real part |k|sqrt(-(D²+4AC))/2, unbounded as |k| grows in the frozen problem.

An explicitly saved scan of the original analytic reference finds

    r = .015542372881355931 m
    z = -.08 m
    direction angle = 1.5882496193148399 rad, from radial toward axial
    D² + 4AC = -74090602.30481993 m²/s²
    growth coefficient = 4303.794904059089 m/s.

Direct matrix eigenvalues independently agree with the discriminant expression.
Uniform density makes A=0 and provides a nonnegative-discriminant control.
The location is interior, whereas the recorded partition failure is at the
axis. This is a counterexample to general hyperbolicity of the isolated
pressure split, **not** a complete explanation of that discrete failure, a
nonlinear well-posedness proof, or a theorem that full warm Hall MHD is illposed.
No unrelated incompressible Hall theorem is substituted for this derivation.

The parent pressure and massless-electron equations follow the species
framework in the [NRL formulary](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp.36–39. The local symbol and counterexample above are this project's
explicit derivation and calculation; they are not quoted from that source.

## Concrete next gate

Derive a **joint** fixed-density Hall-plus-pressure update before treating the
split warm path as generally credible. One candidate uses a common magnetic
midpoint and logarithmic thermal stage: Hall work remains perpendicular to J,
and the paired pressure/thermal work can still cancel globally. Evolve the full
magnetic unknown and U in one nonlinear solve, with consistent pressure and
Hall Poynting/enthalpy face powers and atomic species-positivity rejection.
Do not merely shorten or repeat an ill-behaved pressure substep.

Verify the joint continuum/discrete principal behavior and run actual nonlinear
backgrounds with and without the Hall term. Preserve these exact failed initial
conditions, traces, times and tolerances as controls; add axis-local derivative
and maximum-norm checks. A successful biased-field joint case would not prove
accuracy at an FRC magnetic null or validate electron inertia/kinetic neglect.
The slow warm-wave accuracy and standalone resistive axis defect are additional
requirements. None is resolved by forcing an energy budget or adding arbitrary
heating. Physical collisions, finite mass, vacuum/edge and formation follow.

The joint fixed-density candidate is not a proven cure: Hall dispersion may
degenerate for directions perpendicular to B. If its continuum behavior remains
unacceptable, retain necessary ion/density or finite-mass dynamics from the
parent equations. Do not add a convenient bias or change the physical fixture
merely to remove the observed failure.
