# Experimental nonlinear Hall work pairing

This spatial component is deliberately not enabled in ToroidalMHD: its axis
consistency gate fails. It provides a reproducible discretization and independent
reference with which to repair that failure. The existing working toroidal MHD
solver, its state energy convention, and Stage 2 guard are unchanged.

## State, curl and current

Let x contain corner psi and cell g=Btheta/r. Its magnetic energy H is exactly
the existing cell-poloidal energy plus volume-weighted <r²>g²/(2 mu0). The
linear Faraday map C takes electric coefficients e to xdot:

    psidot = -r Etheta
    gdot = difference_r(Ez)/annular_measure - difference_z(Er/r)/dz.

Etheta is piecewise bilinear in r/z, with zero axis and wall traces. Er/r is
constant in each annulus and linear/periodic in z. Ez is linear in r² in each
annulus, constant in each axial cell, and zero at the outer conducting wall.
Ez at the axis remains free. These are fully filled-cylinder conditions.

Each electric component has a sampling map I onto common tensor-Gauss points,
with physical weights W=2 pi r dr dz and mass M=I.T W I. Component masses are
positive, block diagonal between vector components, and factored once per
geometry. The work-conjugate current is

    M j = -C.T gradient(H),       Jq=I j.

The existing energy is differentiated exactly, including radial centroid
weights, exact <r²>, and folding the duplicated periodic psi endpoint. This is
a variational current, not an independently assumed finite-difference curl.

## Nonlinear Hall projection

At each common quadrature point the full nonlinear Hall field is

    Eq = (Jq cross Bq)/(ne e),
    M e = I.T W Eq,              xdot=C e.

The components of Bq use linear staggered Br/Bz interpolation and Btheta=r g.
Positive carrier density may be supplied at quadrature or explicitly as cell
constants. No density/current floor exists. Then

    gradient(H).T xdot = -j.T M e = -(I j).T W Eq = 0.

This remains true with variable density and finite amplitudes. There is no
Hall Joule heating or additional fluid force. Projected electric samples I e
need not be pointwise perpendicular to Jq; only the raw cross product is.
Both quantities are exposed. The structure-preserving principle is consistent
with the variational Hall cancellation in [Laakmann et al., section 3](https://arxiv.org/pdf/2202.11586).
That paper's incompressible formulation does not prove accuracy of these
specific cylindrical finite-volume energy and electric spaces.

`hall_poynting` separately integrates E cross B/mu0 over faces and takes its
conservative divergence. The local magnetic rate generally differs from this
total-energy flux rate. The difference is recorded, never used to repair
energy or relabeled as physical heating. Radial boundary power is zero; the
periodic axial flux telescopes. No time update or thermal update is installed.

## Independent nonlinear fixture and failed gate

With F=(1-r²/R²)³, the reference is

    psi=B0 r²/2+A r² F cos(kz)
    Btheta=C r F cos(ell kz+.37)
    ne=n0[1+delta F cos(kz+.61)], |delta|<1.

All perturbation field/current traces vanish at the radial wall, and the axis
is regular. Analytic derivatives retain the spatial density gradient inside
curl(E_H). Independent fourth-order differencing checks B, J and E derivatives;
the exact local identity B dot Bdot/mu0+div(E cross B/mu0)=0 holds, including
axis limits. These are prescribed snapshots, not nonlinear time solutions.

The measured physical-volume induction error converges globally, but the first
annulus's gdot coefficient does not: its relative error approaches 21.46%.
The first annulus shrinks in physical volume, so a global norm can converge
while this coefficient remains wrong. Zero-toroidal/poloidal perturbation,
carrier sampling and quadrature controls are retained. Exact global Hall work
does not rescue this failure. See HALL_WORK_VERIFICATION for the measured data.

The next change must reconcile regular-axis polynomial moments, current and
electric sampling, and the magnetic energy convention. Any change to energy
must be explicit in pressure recovery and Poynting accounting; a fitted axis
multiplier or post-step energy reset is not an accepted solution. After that,
verify local transfer and advance all compressible wave branches with actual
time/mesh refinements. The current SSPRK2 fluid step has |R(i x)|²=1+x⁴/4 for
an undamped imaginary mode, so adding a nondissipative Hall RHS alone is not a
stable Hall time integrator. A suitable time scheme, Hall scale limit or
converged implicit solve is required and must be measured independently.
