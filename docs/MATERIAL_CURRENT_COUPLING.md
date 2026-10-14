# Finite-mass current on a moving material boundary (2026-10-14)

This advances the ideal magnetic component from c4c32463 to independent evolving
current/field coordinates, simultaneous with the moving fluid. The implemented
axisymmetric invariant sector has toroidal B and meridional bulk/carrier flow.
Both species are isentropic with equal constant entropy and eta=0. It retains
physical electron/ion masses, Hall/self-drift, actual changing-density inertia,
fluid vorticity/stretching, magnetic force and current stress. It is not yet the
full thermal, poloidal-field, vacuum/coil and FRC formation system.

## Physical equations and boundary

Let M=mi+me, mu=mi*me/M, alpha=me/M, beta=mi/M, delta=beta-alpha, and
w=vi-ve=J/(en). The parent bulk/current equations in this sector are

    D_t n = -n div(v),
    M n D_t v = J x B - grad(pi+pe) - div(mu n w w),
    E+v x B = delta w x B + [alpha grad(pi)-beta grad(pe)]/(en)
              + (mu/e)[D_t w+(w.grad)v-delta(w.grad)w].

The field/current energy includes integral B²/(2mu0)+mu*n*|w|²/2. The fully
material covariant drift basis gives

    E* = delta w x B* - (mu/e) w x curl(v) + (mu/e) G wdot + grad(Phi),
    B* = B + (mu/e) curl(w),  E*=E+v x B.

The barotropic pressure and self-drift scalar potentials have zero compatible
curl. They are excluded from the returned reduced electric field, not treated
as absent physical pressure. Both species pressures drive actual fluid motion.
Constant species entropy makes their relative entropy transport identically zero;
equal temperatures make thermal exchange zero. Eta=0 is an explicit mathematical
limit, not a Coulomb-transport recommendation at the chosen n,T.

The modified frozen flux associated with carrier inertia is discussed by
[Lingam, Morrison and Tassi](https://arxiv.org/abs/1410.8152), and Hall/inertial
extended-MHD relations by [Lingam, Morrison and Miloshevich](https://arxiv.org/abs/1504.08293).
The dimensional finite-mass coefficients and expansion reference below are
derived here from the repository's two-fluid parent equations; these papers do
not validate this mesh, dilute boundary or laboratory regime.

Both cylinder ends and its side move with the fluid. Reference Btheta has zero
traces on all three surfaces, preserved by the physical Piola map. Regular
axisymmetric toroidal magnetic vacuum has Btheta=0, so the exterior has exactly
zero field, matter and conductivity in this sector. No periodic axial wrap,
sealed cap, numerical atmosphere or density floor is present. Penetrating
poloidal magnetic fields and actual external coils are still missing here.

## Moving compatible discretization

The existing nodal/histopolant polynomials in r0²/2,z0 give covariant electric
basis G and contravariant magnetic basis B. The raw incidence annihilates nodal
gradients. A reference trace null space imposes the zero toroidal field. Actual
physical quadrature metrics are assembled at every stage. In particular, the
old fixed-reference Faraday projection is not used after deformation.

With fixed material particle weights dN,

    Mn = integral dN G^T G,
    H = integral Bbasis^T Bbasis dV / mu0,
    D = -integral Bbasis^T curl(G) dV / mu0,
    e Mn w = -D^T b,  H bdot = D E*.

J=en Gw therefore disappears with carriers; it is not formed by dividing a
cellwise curl by a floor. Curl(G) includes the complete reference-compatible
map. Weak Ampere includes field jumps. The compatible gradient annihilates its
current charge work, but pointwise charge continuity and local face transport
are not thereby proven.

The actual finite-inertia solve is

    (H + mu/e² D Mn^-1 D^T) bdot = D Mn^-1 rhs,
    rhs = f - mu/e² Ddot^T b - mu/e Mndot w,
    f = integral dN G^T [delta w x B* - (mu/e) w x curl(v)].

Both geometry terms use the evolving shape and velocity. Magnetic force is the
fixed-flux energy shape derivative; bulk current stress is
integral mu dN ww:grad(delta x). Their geometric work cancels that in Wb+Kj.
The same RK4 step advances material coordinates, momenta and independent magnetic
coefficients. No time step changes, diffusion, energy repair or altered masses
are hidden in this model.

## Independent expansion reference

For a prescribed homothety lambda=1+t/tau, n=n0/lambda³ and zero vorticity,
linearized extended MHD transports B* flux. The finite-cylinder mode

    Btheta0=B0 J1(kr r) cos(kz z),
    kr R=first zero of J1, kz=pi/(2L),
    chi=mu*(kr²+kz²)/(mu0 e² n0)

has independent material amplitude

    b(t)/b(0)=(1+chi)/(1+chi lambda),
    Bphysical/B0 = [b(t)/b(0)]/lambda²,
    Jphysical/J0 = [b(t)/b(0)]/lambda³.

The runs choose chi=1 with physical Ar+/electron masses and retain nonlinear
Hall/self-drift terms at B0=1e-12T. Their characteristic nonlinear parameter
e*B0*tau/mu is about1.76e-7, below measured spatial errors; this is a weak-field
reference, not a claim that Hall is removed from the solver. At lambda=2 the
physical amplitude is1/6, whereas ideal flux freezing gives1/4.

Current energy is comparable to magnetic energy. For this mode,
Wb=Wb0*[b/b0]²/lambda and Kj=chi*Wb0*[b/b0]². The geometric power satisfies
(Wb+Kj)'=-(lambda_dot/lambda)*(Wb+2Kj). Recorded work uses trapezoidal integration,
so its time-quadrature error is separate from the RK4 solution error.
Prescribed fluid motion does external work: its reported total-fluid energy
change and unforced virial residual are not conservation checks for that mode.

## Actual freely moving model

The free run starts with
n0=nc*[(1-r²/R²)(1-z²/L²)]^(3/2) and Te0=Ti0=[(1-r²/R²)(1-z²/L²)]eV.
Each material entropy is constant; density and pressure vanish at the side/ends.
The compact initial Btheta=(2e-4T)*(r/R)*[(1-r²/R²)(1-z²/L²)]³ is projected
onto the actual boundary-admissible magnetic space. Projection error is part of
the initial data error, not a change to the continuing field equation.
R=35mm,L=70mm,nc=1e18/m³; duration R/sqrt(10e/M)=7.12177us.

For the unforced full finite-mass sector, the independent scalar virial is
I''=4Kbulk+4Kj+4U+2Wb=4E-2Wb. Thus the full current kinetic energy belongs in E.
Its integrated check uses the actual evolving Wb and a time quadrature; no
conservation discrepancy is fed back to the trajectory. A cold dilute edge and
eta=0 do not establish laboratory validity. General transport and ionization
remain required regardless of numerical conservation.

## Linear reconstruction failure and repair

The original fine free run stopped at attempted step554, preserving553 accepted
states. Its Faraday backward residual1.031048e-10 exceeded the unchanged1e-10
limit while Ohm was2.36e-16. This is cancellation in reconstructing E after
eliminating the mixed equations, not evidence that a looser tolerance is valid.

The revised solver refines the original two rows
H*bdot-D*E=0 and (mu/e²)D^T*bdot+Mn*E=rhs using the same Schur operator, up to
four correction passes. These corrections solve the linear equations; they do
not modify accepted physical states, forces or energy. The original Faraday/Ohm
checks remain unchanged. The failed source/output remains immutable. The repeated
fine trajectory completes all 2000 steps with maximum Faraday residual 3.04e-13
and Ohm residual 7.76e-16. This resolves the observed linear reconstruction failure;
it does not resolve spatial current error.

## Initial field and carrier representation

Magnetic-only initial projection gave current kinetic energies 84.73% and 29.21%
above the independent continuum integral on 4x8 and 8x16 cells. The corresponding
magnetic energy errors were only -0.860% and -0.0493%. These original complete
trajectories are preserved; a good magnetic energy alone hid a poor current.

The new `magnetic_current` initial projection minimizes

    0.5 integral |Bbasis b-Bexact|² dV/mu0
      + 0.5 mu integral |G w(b)-wexact|² dN,
    e Mn w(b)=-D^T b.

The mutually consistent target B and J are independently differentiated compact
formulas. The same positive matrix H+mu/e² D Mn^-1 D^T defines this initial
minimization. It never adjusts a continuing state or fits an energy total.
An independent 80x80 continuum integral gives Wb=5.22144522e-8 J and
Kj=4.80416085e-8 J. The tradeoff is measured, not hidden:

| Initial error, new projection | 4x8 | 8x16 |
|---|---:|---:|
| Magnetic field, relative L2 | 16.74% | 5.212% |
| Current, relative L2 | 32.23% | 16.15% |
| Drift, particle-weighted relative L2 | 32.08% | 25.37% |
| Magnetic energy, signed relative | -15.48% | -6.055% |
| Current kinetic energy, signed relative | +3.490% | -0.1517% |
| Acceleration RMS / (C/R), C=10e/M | 0.002869 | 0.0008402 |

The independent acceleration includes pressure, JxB and the divergence of the
full drift stress, with cylindrical connection terms. Initial B errors increase
relative to magnetic-only projection (9.274%/2.221%), while current and drift
errors improve. Near-correct integrated Kj does not imply an accurate local J.

## Completed trajectories and measured limits

The two prescribed-expansion runs complete 400 steps over 1 us. Final magnetic
RMS errors normalized by the initial field scale are 6.5272e-4 and 1.3378e-4;
current errors normalized by the initial current scale are 1.0013e-3 and
3.7145e-4. Mesh doubling improves them 4.88-fold and 2.70-fold. Both follow the
finite-inertia amplitude, which is distinct from ideal flux freezing. Maximum
Faraday residual is 1.79e-14. These are mathematical comparisons, not experiments.

The new-projection free runs both complete 2000 steps over 7.121769 us:

| Full free trajectory | 4x8 | 8x16 |
|---|---:|---:|
| Maximum relative total-energy change | 8.57e-14 | 6.44e-10 |
| Maximum independent virial error | 3.55e-12 | 3.21e-11 |
| Final magnetic energy, nJ | 28.1488 | 15.4129 |
| Final current kinetic energy, nJ | 52.3884 | 63.9767 |
| Final peak magnetic field, uT | 33.3805 | 26.7792 |
| Final peak current, A/m² | 3153.06 | 14280.03 |
| Final peak relative drift, m/s | 310655 | 419638 |
| Wall time, seconds | 121.1 | 432.3 |

Maximum Faraday/Ohm residuals are 2.63e-13/1.25e-15; compatible weak charge work
is below 6.78e-17. **The free trajectory is not spatially converged:** final Wb
changes about 45%, Kj about 22%, and peak J about 4.5-fold. Conservation and the
virial identity do not establish accurate local field/current structure. The
final spatial plot visibly retains unresolved structure; it must not be used
as a laboratory prediction or a converged formation result.

The existing conditional Coulomb helper samples maximum density, maximum drift
and minimum density in initial/final states. For the finer run it gives roughly
34–190 electron collision times over the simulated duration, flagging eta=0 as
a formal limit. Final drift/electron-thermal-speed ratios reach 12.2 at the
maximum-drift sample and 24.0 at the minimum-density sample. These estimates
assume the helper's Maxwellian/weak-coupling regime and use body radius, not a
resolved gradient length. They do not certify a transport closure at large
drift. Do not silently install Spitzer resistivity as an edge-model repair.

Pointwise current continuity, normal relative-drift boundary traces, general
species entropy/heat transport and kinetic/ionization edge behavior remain
unresolved. All exterior matter and conductivity remain exactly zero. This
toroidal-field fixture has no FRC poloidal topology or external coil drive.

## Evidence and next physical work

The inspected [scientific figure](../outputs/material_current_comparison/material_current_expansion.png)
and its PDF show the independent expansion law, field/current errors, free-run
energies and final field/current maps. The adjacent `comparison.json` includes
initial representation errors, historical comparisons, conditional validity
samples and source-summary hashes. Accepted states and frozen sources are in
`outputs/material_current_reference_r{4,8}_t400` and the five
`outputs/material_current_free_*` directories. The original failed fine run is
retained alongside the completed mixed-row-refined and new-projection runs.

Scientific sessions 15237, 1742, 50963 and batch 62098 are terminal with exits 0
collected. Original failure 94090 is terminal 1; observers 95334 and 43866 are
terminal 0. No worker is live. No unit tests or new verification framework were
added or run. The older interactive 3D viewer is unchanged.

Next improve the actual spatial field/current representation using the existing
polynomial degree and quadrature controls, then compare a complete trajectory
with a necessary time-step control. Do not infer convergence from energy totals.
Follow with general species entropy/heating and poloidal-field/moving-vacuum/
finite-coil coupling. Preserve the fuller fixed-footprint finite-coil solver;
this restricted material sector does not replace it. `frc-vacuum.3`, its parent
and the full FRC goal remain open.
