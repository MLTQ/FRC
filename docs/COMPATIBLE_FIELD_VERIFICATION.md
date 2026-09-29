# Compatible magnetic/current verification and harmonic precision

The new higher-order field component resolves the previously diagnosed slow
linear magnetic-rate error under mesh refinement. **343 tests pass in 118.34 s**.
All 25 affine-field cases and 25 preserved rounded-total comparisons pass the
independent source/field/equation/physical audit. This verifies instantaneous linear components, **not** the complete
nonlinear moving thermal/energy update, a formation trajectory or confinement.

Equations, basis, exact harmonic split, boundaries, absent physical/numerical
diffusion and required coupling work are in
[COMPATIBLE_FIELD_NUMERICS](COMPATIBLE_FIELD_NUMERICS.md).
The earlier nonlinear smooth solver is unchanged; its pressure improvements
and remaining full coupled errors in SMOOTH_FINITE_MASS_VERIFICATION still apply.

![Compatible field convergence and retained precision failure](../outputs/compatible_wave_verification/accuracy.png)

The standalone scientific figure was visually inspected; its final hash is
retained with the report.

## Independent checks and declared domain

Five component tests check nodal polynomial values/derivatives, true interval
integrals, nonuniform/periodic commuting derivatives, independent cylindrical
polynomial B and curl(B)/mu0, physical magnetic energy, work, original Faraday
signs/metrics, zero curl of a general potential, conservative variable-density
velocity reconstruction, and a separately assembled original three-block
finite-inertia solve. The polynomial current respects the actual mixed wall
trace; an incompatible nonzero wall Jtheta is not quietly fitted away.

Two further tests check exact zero current of a uniform physical bias and
polynomial currents at amplitudes 1,1e-4,1e-8,1e-12 beside it. Perturbation energy
is normalized by its own squared amplitude. Independent physical local energy
retains the nonzero spatial cross term, which integrates to zero only globally.
The initial random-work assertion used an unsuitable near-cancelled net
denominator; its retained failed log is replaced by an explicitly bounded
noncancelled floating-point work check. Separate local physical accuracy tests
are unchanged. No physical accuracy tolerance was relaxed to pass that identity.

The exact finite-mass compressible Bessel/Fourier oracle is independently derived
from the original linear species equations. Main cases use R=.035 m, L=.24 m,
n0=1e18 m^-3, B0=.005 T, physical Ar+/electron masses, reference Te=5 eV,
Ti=1 eV, gamma_e=5/3 and gamma_i=1.4. Equal-energy wave amplitude is 1 m/s.
The three angular frequencies are approximately 2807.94, 545649 and
53924707 s^-1. These are linear perturbations about a constant background;
the bulk momentum input is the analytic cell average of rho0 v. Both J and
v used in the numerical force are reconstructed, never substituted from the
analytic reference. The reference supplies only initial data and error targets.

The protocol covers degree5, nr8/16/32/64, nz=2 nr for all three branches;
four degree3 slow controls; four explicitly matched perturbation-only .01-amplitude
linear controls; quadrature orders 8/10/12 at nr32; a hundredfold weaker full-bias
wave; original n0=1e20; and opposite slow propagation. Main states retain the
full physical bias as a separate harmonic. Stored increments, full rounded
export, velocity, forcing, potential, current and all rate/electric arrays are
preserved. The exact Bessel root gives zero wall increment explicitly.

## Magnetic-rate convergence

Errors are independently integrated physical relative volume L2, divided by
the actual magnetic-rate signal. They are not normalized to the static B0.

| nr | Slow | Intermediate | Fast |
|---:|---:|---:|---:|
| 8 | .789352 | .00143441 | .0000862078 |
| 16 | .0141969 | .0000270477 | .00000268075 |
| 32 | .000190057 | .000000423949 | .0000000820072 |
| 64 | .00000236657 | .00000000667951 | .00000000250544 |

Slow nr64 error is about **0.000237%**, with dimensional RMS 3.38e-12 T/s,
sampled maximum 2.21e-11 T/s and first-annulus relative L2 8.55e-7.
The intermediate first-annulus rate error is 1.97e-7, larger than its whole-volume
6.68e-9. Retaining axis/max metrics matters even when global errors are small.
These observed rates do not establish a universal order for nonlinear coupling.

For the exactly matched weak, perturbation-only linear inputs, the previous
representation gives slow rate errors 1004.98/263.069/66.8838/16.8431.
The compatible degree5 values are .789352/.0141969/.000190057/2.36657e-6.
This is a comparison of **linear magnetic components**; it must not be used to
claim the earlier complete nonlinear trajectories have already been corrected.

Degree3 slow errors are 31.3134/1.75203/.0922928/.00502551. Degree5 materially
improves this cancellation at the same node count. No fitted current, modified
electron mass, changed density, shorter time window or exact forcing field is
responsible for the improvement.

## Separate B, J, E and current-rate accuracy

Slow-mode errors in the initial physical B, Ampere J, full E and J_t:

| nr | Initial B | Initial J | Full E | J_t |
|---:|---:|---:|---:|---:|
| 8 | 9.65402e-5 | .000438729 | .000553480 | 3.86867 |
| 16 | 3.29637e-6 | .0000272033 | .0000329164 | .129285 |
| 32 | 9.78096e-8 | .00000137804 | .00000163204 | .00330316 |
| 64 | 2.81784e-9 | 6.46807e-8 | 7.58444e-8 | .0000773307 |

Current-rate accuracy is less favorable than magnetic-rate accuracy; the coarse
grids do not resolve it adequately. Finest slow J_t RMS/max errors are
.00990395/.115208 A m^-2 s^-1, with first-annulus relative error 8.04e-7.
The finest intermediate/fast J_t errors are 9.17e-8/4.57e-9. Initial-current
errors for those branches are 4.67e-8/6.47e-8. All three retain full dimensional
and sampled maximum norms in the records.

Slow velocity-reconstruction errors are 6.49e-6/1.33e-7/2.01e-9/2.61e-11.
Velocity uses actual conservative annular averages, independently checked with
variable-density polynomial data. This still does not supply the full moving
vorticity, local kinetic work or conservative thermal coupling.

## Bias roundoff and explicit affine storage

The first high-order implementation retained one rounded total magnetic array.
With full B0 its slow nr8/16/32/64 rate errors were .789419/.0159158/
.0117924/.0245827. Refinement eventually made the error worse. A hundredfold
weaker nr32 wave gave **1.17396**, while the perturbation-only control gave
.000190057. This initial study remains in `compatible_wave_rounded_attempt`;
its fields, source, failed accuracy and complete records are not overwritten.

At nr32, evaluating current after subtracting the known harmonic from the same
already-rounded input reduced the exploratory magnetic-rate error to .000665870.
Initializing/storing the increment directly gave .000190057. This distinguishes
operator cancellation from information already lost in total-field storage;
subtraction cannot recover missing digits.

The final explicit affine field stores the physical bias and increment separately.
Its full-bias and hundredfold weaker nr32 errors are .000190056627 and
.000190056632. The finest full-bias and matched weak controls likewise agree
within measured roundoff. Uniform bias retains its field, energy, force and
local cross term; only its mathematically exact zero-current contribution is
removed from floating-point differentiation. This is an exact representation
change for the stated boundaries, not an imposed equilibrium correction.

The lossy rounded export still has nonzero original-Ampere error when replayed
relative to its tiny current: for slow nr64 this is 4.72e-6. Its un-cancelled
original backward error remains separately recorded. The affine representation
does not claim that export preserves every increment or that arbitrary static
fields can be treated as a uniform harmonic.

## Controls, equations and physical energy

At nr32, assembly quadrature orders 8/10/12 give slow magnetic-rate errors
.000190056688/.000190056627/.000190056691. Opposite propagation gives
.000190056648. Original-density n0=1e20 gives 9.67998e-5, with current-rate
error .00193148. The harder original physical regime has not been replaced
by this denser control.

Independent physical norms use root-product polynomials and orders 10/14,
without production basis solves, sampling or energy matrices. Both full and
small-increment magnetic energy, physical current energy, projected bulk/Hall
force, pressure potential and original Ohm/Ampere/Faraday equations are audited.
The largest independent physical energy discrepancy, including normalization
by the small increment energy itself, is 2.75e-14 relative. Norm quadrature
changes are at most 7.11e-14 absolute in relative L2 units. Original equation
backward errors reach 1.83e-12; separate full-field Ampere backward error is at
most 1.31e-16. The rounded export's worst current-relative Ampere replay error
is .000129808, so those backward and small-signal measures are not interchangeable.
The unchanged small-signal solver's maximum residual is 7.75e-11 at tolerance
1e-10. After restoring the full pressure potential, full-array Faraday error
relative to its tiny rate reaches 5.88e-10; the original backward error remains
distinct. No claim is made that this is below the small-signal solver tolerance.

No density, momentum or temperature is evolved by these rate cases. There is
therefore **no full-trajectory energy or entropy claim** in this study, and no
inference that small linear residuals prove physical accuracy.

## Reproducibility and required next action

Rounded-attempt source:
`485edaa1fdb9bffcb76b2a4e9196e3e4b1ddc071b20f319b4c3062b42ebae667`.
Affine measurement source:
`8cfb1a2f89f7ee58fadae2e47ae9e4fbb5c49037face19da117518188ed8d720`.
Both retain source snapshots and every case's exact configuration, input/rate
arrays, hash, norm and conditioning result. A separate report/helper/test archive
records the independent audit. All old smooth finite-mass sources remain unchanged.
The bias-cancellation probe has its own input, output and helper-source archive.
`checkpoint.json` records final source/report/figure/test/document hashes and
terminal job statuses. No remote repository or publication is configured.

Run fresh measurements with
`python -m examples.compatible_wave_verification --directory <new-path>`.
The driver rejects output-directory reuse. Rebuild the read-only comparison
with `python -m examples.compatible_wave_report`. The report, scientific figure,
test/worker logs and final hashes belong to `outputs/compatible_wave_verification`.

Next derive the variable-density physical current-energy/pressure split and
complete common work/particle transport on this field representation. Retain
the actual mass-flux n_t and local harmonic cross energy. Verify the parent
nonlinear manufactured equations and full contact/wave trajectories before
adopting it in SmoothFiniteMassMHD. The heated free edge/exact vacuum, physical
transport/ionization, shock/entropy behavior, real ends/coils, converged FRC
formation and interactive experiments remain unfinished. Axisymmetry alone
would establish neither 3-D stability nor particle confinement.
