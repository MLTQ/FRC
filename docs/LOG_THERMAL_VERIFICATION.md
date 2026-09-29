# Physical thermal metric: measured accuracy and limits

The new physical log-entropy component corrects the previous nonconverging local
axis rate in matched smooth tests. All **18 rate cases and nine full-transit
transport cases** completed and passed source, array, original-equation and
independent physical audits. The complete suite passes **361 tests in 121.42 s**.
This is component verification, not a new full moving-plasma or formation result.

The equations, boundary assumptions and density/energy contract are in
[LOG_THERMAL_NUMERICS.md](LOG_THERMAL_NUMERICS.md). Boundaries are periodic z,
regular axis, fixed poloidal wall flux and natural Btheta(R)=0. The tests use
strictly positive prescribed conservative density; no vacuum atmosphere, floor,
fitted current/mass or energy reset is introduced. The spatial Galerkin component
has no added viscosity or limiter. RK4 truncation and unresolved spatial
transport produce the separately measured energy/proxy drift below.

## Protocol and independent references

The rate protocol uses degree 3/5 at nr=8/16/32/64, nz=2nr, quadrature 8/10/12,
physical density 1e16/1e18/1e20 m^-3, opposite current, constant-entropy and
uniform-temperature controls. Both species use their physical finite masses;
electron gamma=5/3 and ion gamma=1.4 exercise separate closures. Actual
conservative density and the supplied continuity rate are unchanged from the
variable-current initializer. Analytic entropy is not fitted to the discrete
current. The current-only pressure pair and changing-density induction are
combined with their declared exchange with the missing bulk dynamics.

Independent root-product scalar/field bases and primitive-polynomial density
reconstruction evaluate physical volume norms and local energies at quadrature
12/16; the audit repeats them at order 14. Polynomial pressure and an independently
assembled scalar mass solve test the metric. Directional energy differentiation
checks both entropy and actual conservative density. The transport reference
uses a complex-step PDE check, closed Bessel/radial integrals and an independent
matrix exponential for fourth-order time convergence. These checks complement
the algebraic work pairing; they are not inferred from it.

## Matched correction of the local axis rate

All entries below are relative L2 errors. “Axis” here means the first four
nodal radial rows, each normalized by the exact signal in those rows. The old
log rate is its proxy rate divided by its nodal proxy; the new rate is the
primary log-state rate. Both compare with the same analytic nodal target.

| nr | Old electron axis log rate | New electron axis log rate | New physical proxy rate | Old physical proxy rate |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0.309416 | 1.02400e-3 | 2.37926e-3 | 0.746408 |
| 16 | 0.374235 | 2.96189e-5 | 4.76038e-5 | 0.411293 |
| 32 | 0.387771 | 5.42997e-7 | 8.16064e-7 | 0.229013 |
| 64 | 0.390908 | 8.84960e-9 | 1.30658e-8 | 0.117630 |

The continuous proxy comparison samples both representations in the same
physical nodal basis and volume norm. It is distinct from the previous report's
6.35% whole-domain lumped nodal proxy-rate error. Likewise the four nodal rows
are distinct from the physical first-annulus norm used for transport below.
Shrinking the axis region is not hiding a fixed local error: its locally
normalized new rate converges, while the old one does not.

| nr | Electron physical log rate | Electron pressure E | Ion pressure E |
| ---: | ---: | ---: | ---: |
| 8 | 1.84763e-3 | 2.50779e-3 | 2.48762e-3 |
| 16 | 3.75480e-5 | 7.65851e-5 | 6.47152e-5 |
| 32 | 6.43018e-7 | 1.77143e-6 | 1.61712e-6 |
| 64 | 1.02909e-8 | 4.37041e-8 | 4.21329e-8 |

Degree 3 retains lower convergence: finest electron proxy/E errors are
1.55761e-6/8.29082e-6. Quadrature and density controls retain the degree-5
improvement; current reversal leaves those error norms unchanged. Full absolute,
sampled maximum, physical axis, local energy, pressure and temperature norms
remain in the raw report rather than being replaced by these selected L2 values.

## Work, entropy and full electric-field precision

Across the rate protocol the maximum pressure-work backward defect is
1.85e-14, physical entropy-rate defect 9.68e-15 and declared total-energy exchange
defect 1.37e-14. Independent local energy/rate evaluation differs by at most
6.83e-15. Scalar load residuals are at most 1.99e-13; induction residuals are at
most 1.91e-11 at their unchanged 1e-10 tolerance. Original full-equation backward
errors are below 3.16e-13. These normalized cancellation/residual measures are
not continuum errors.

The physical log entropy is the conserved quantity. The exponential proxy
integral is not another exact discrete invariant. Its signed rate defect divided
by the integral of the absolute proxy rate decreases from 2.62e-8 to 3.58e-12
and then floating-point scale in the default refinement. No shock entropy
inequality or exact finite-step energy conservation follows from these results.

Recombining the removed gradient with the reduced electric field loses some
curl precision. In the finest uniform-temperature case, replaying Faraday from
the full rounded E differs from the stored rate by **0.00242218 in the physical
magnetic-energy norm**, or **6.15e-10 T/s RMS**. The separate raw coefficient-space
relative difference is **0.0203616**; the two norms must not be interchanged.
The largest physical RMS discrepancy over the protocol is 2.71e-6 T/s in the
n=1e16 case, where its rate-relative value is 3.05e-10. Thus small backward error
does not assert uniformly accurate weak-signal reconstruction from exported E.

## Constant entropy and uniform-temperature controls

Constant nodal entropy gives exactly zero stored entropy rate, reduced pressure
field and joint magnetic rate. Its restored full gradient has a roundoff curl
of 1.93e-12 T/s RMS; a relative rate error is undefined for the zero target.
Pressure E itself still has spatial approximation error.

At varying density, uniform temperature is generally not exactly representable
in the finite log-entropy polynomial space. The continuum pressure battery is
zero, but the computed stored rate is small and nonzero:

| nr | Initial electron T relative error | False B_t RMS (T/s) | Sampled maximum (T/s) | First-annulus RMS (T/s) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 6.85086e-5 | 4.16574e-2 | 1.45512e-1 | 8.58778e-3 |
| 16 | 1.41499e-6 | 9.21031e-4 | 3.33532e-3 | 9.38127e-5 |
| 32 | 2.39292e-8 | 1.58422e-5 | 5.85793e-5 | 8.03015e-7 |
| 64 | 3.81638e-10 | 2.53764e-7 | 9.43100e-7 | 6.42308e-9 |

This is converging false response, not exact preservation of an isothermal state.
For general nonbarotropic rate cases no analytic joint B_t oracle is supplied;
their absolute rate is not mislabeled as magnetic error.

## Complete prescribed-current transport

The independent wave has n=n0(1+0.3 r²/R²), Btheta=r G(1-r²/R²)² and
s=s0+A cos(k[z+Jz t/(e n)]). Its sign-changing axial current has zero net flux;
all radial wave phases complete the same elapsed axis-transit time
**1.208012722e-6 s**. Density and current are prescribed throughout. This is
neither a force-balanced plasma nor a solution of the full coupled equations.

Nine runs include nr=8/16/32, half dt, higher quadrature, reversed current,
larger amplitude, gamma=1.4 and degree 3. Errors use the entropy/temperature wave
signal as denominator, not the large background entropy.

| nr | Full-transit entropy error | First-annulus entropy error | Temperature-wave error | Signed internal-energy drift | Signed proxy drift |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 8.46755e-2 | 1.64112e-2 | 8.61826e-2 | 9.56856e-6 | 1.78566e-6 |
| 16 | 2.03068e-3 | 2.23039e-5 | 2.12812e-3 | -3.59613e-9 | 9.83506e-10 |
| 32 | 9.31941e-5 | 9.65239e-9 | 1.02875e-4 | -1.46738e-12 | 3.60756e-12 |

These drifts are relative to each discrete initial state. Initial physical
energy/proxy errors against the exact integrals are separately -5.02e-7,
-8.00e-9 and -1.26e-10. At nr=16, halving dt changes the final entropy field
by 1.73e-8 of its wave signal, far less than the spatial error. Quadrature and
current reversal agree to near floating-point precision. The amplitude-1
control retains the same relative entropy error but has energy/proxy drift
6.26e-8/1.06e-7; nonlinear observables are not silently assumed amplitude-free.
Degree-3 entropy error is 0.00700418 at nr=16.

Maximum physical entropy drift is 2.22e-16 across all nine runs. Independent
snapshot integrals agree within 5.23e-15. All accepted pressures are finite and
positive, with atomic rejection tests for unrepresentable trial states. The audit
replays initialization and the final accepted RK4 step bit for bit, and independently
recomputes physical norms/integrals for all five saved snapshots. It does not
claim bitwise replay of every earlier step.

## Reproducibility and next action

Final source identity:
`b0dba42d2dd23742246776328413c3914ecfebaadf158e72d9b157c10d6485e7`.
Records, raw arrays, source/companion archives and audits are in
`outputs/log_thermal_verification` and `outputs/log_entropy_evolution`.
The report helper has a separate archive/hash outside the frozen calculation
closure. Test and worker logs are retained; all workers are terminal. The
[six-panel figure](../outputs/log_thermal_verification/verification.png) was
visually inspected. The prior variable-current, compatible-wave and nonlinear
smooth source identities remain unchanged.

Next derive a common full moving state and energy contract. Two nodal log-entropy
fields cannot simply replace the earlier conservative cell-energy/ion-pressure
partition. Continuity, bulk entropy advection, conservative velocity/vorticity,
momentum/drift stress and all pressure/current/enthalpy/Poynting work must use
compatible physical reconstructions. Verify independent parent-species rates
and full wave/contact trajectories before adopting the new representation.
Heated free edges and exact vacuum, physical transport/ionization, real end/coil
coupling, shocks and kinetic validity, converged formation and the UI remain
open. Axisymmetric accuracy establishes neither 3-D stability nor confinement.
