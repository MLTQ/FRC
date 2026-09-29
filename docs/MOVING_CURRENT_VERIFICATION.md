# Moving finite-mass induction: measured evidence

The full suite passes **322 tests in 107.36 s**. All **19 instantaneous rate
cases** complete and pass the independent archive/physical audit. This verifies
the moving-current equation and a continuum momentum/energy contract. It is
not yet a conservative moving finite-mass time integrator or a formation run.
Equations, boundaries and numerical approximations are stated in
[MOVING_CURRENT_NUMERICS](MOVING_CURRENT_NUMERICS.md).

## Independent parent species

The manufactured cylinder has R=0.035 m, L=0.24 m, physical Ar+/electron
masses, n0=1e18 m^-3, Te=5 eV, Ti=1 eV, gamma_e=5/3 and gamma_i=1.4.
Nonuniform n, both pressures, all three bulk velocities, poloidal/toroidal B
and time-dependent compression are analytic. Exact particle fluxes satisfy
continuity and analytic E satisfies Faraday. Each species' prescribed force
and heat are computed from its original equations, before reduction; no
forcing is fitted to the grid. Original-density n0=1e20, equal physical masses
and consistent massless cases are separate verification controls. These
inputs are not a claim of physical regime validity.

Continuity, bulk momentum, Ohm reduction, the B* transformation, Faraday and
direct species energy-flux sums agree to at most 5.09e-15 relative, including
samples approaching r/R=1e-6. Independent fourth-order physical differentiation
of full local energy flux gives residuals 4.6518e-10 → 2.9062e-11 → 1.8624e-12
for physical argon as the spacing halves. The equal-mass control gives
3.8637e-9 → 2.4149e-10 → 1.5281e-11.

Omitting mixed bulk/current energy transport leaves a nonvanishing local
defect: 9.68227e-4 for argon and 6.19090e-3 for equal masses, relative to the
manufactured power norm. These are not total-energy drifts. In the argon
fixture, drift-stress divergence reaches 20.8468 N/m³; changing-density and
velocity-curl electric terms reach 2.11425 and 0.656895 V/m. At n0=1e20 these
contributions are 100 times smaller; consistent massless controls set them
exactly to zero.

## Spatial rate refinement

Physical volume-relative L2 errors compare with analytic B_t, J_t and full E,
using independent interpolation/quadrature. nz=2 nr. The denominator is the
actual rate/field, not the stationary magnetic bias.

| nr | B_t error | J_t error | E error | First-annulus B_t error |
|---:|---:|---:|---:|---:|
| 8 | 0.217882 | 0.275820 | 1.54929 | 0.0744755 |
| 16 | 0.0968733 | 0.0718900 | 0.287859 | 0.0182798 |
| 32 | 0.0469688 | 0.0241030 | 0.0589906 | 0.00454100 |
| 64 | 0.0232974 | 0.00985879 | 0.0161778 | 0.00113331 |
| 128 | 0.0116252 | 0.00459108 | 0.00643579 | 0.000283210 |

Finest first-annulus J_t/E errors are 0.00711149/0.00708714. B_t approaches
first order in volume and second order in the first annulus; J_t/E remain
approximately first order in the first annulus. A small equation residual
does not remove physical interpolation/truncation error. At nr128, absolute
RMS errors are 29.6501 T/s, 1.26238e9 A/(m² s) and 0.167987 V/m; sampled
maxima are 102.833 T/s, 8.03060e9 A/(m² s) and 0.989222 V/m. Maxima depend
on sampling; order-8 and order-12 values are both archived.

The n0=1e20 control reaches B_t/J_t/E errors
0.0232367/0.00920142/0.0119987 at nr64. The n0=1e18 massless comparison is
substantially less accurate in J_t: its nr64 volume error is still **0.403719**,
and its latest observed order is only 0.588. This limitation remains visible;
passing the solve is not adequate current-rate accuracy. It is not evidence
of continuum ill-posedness or proof of an earlier split failure's cause.

![Physical rate refinement](../outputs/moving_current_verification/refinement.png)

The figure was visually inspected. Raising norm quadrature 8→12 changes the
relative volume/axis norms by at most 1.33e-15 absolute. Production quadrature
5→7 changes them by at most 4.41e-10; tighter linear tolerance by at most
1.68e-13. Two extra compression phases, an equal-mass case and a Galilean-
translated analytic finite-inertia wave pass. The translated wave retains
its explicit bulk constraint; it is not an unforced compressible-wave test.

## Moving terms, work and preserved failures

At nr32, separately dropping the density, vorticity or bulk B* term changes
the magnetic-rate energy norm by 1.89967%, 0.146032% or 28.0284%. These
deliberately incomplete equations are diagnostic controls, not alternate
physical models. Their fields and equation replays are retained.

Every original equation residual meets its requested tolerance: 1e-10 in the
main cases, 1e-12/2e-12 in nr16/nr64 controls. The largest recorded residual
is 1.5898e-11. Independent full-E Faraday replay is at most 1.3718e-11.
Changing-density magnetic/current energy derivatives balance forcing and the
explicit derivative of 1/n to at most 1.3690e-10 of component power scales.
Both differentiation intervals are retained. These are tangent checks, not
finite-time conservation. The finest density-energy contribution is -0.879098 W.

The initial protocol stopped at target nr32 when subtraction of nearly equal
global energies lost precision. Its nr128 solve also rejected residual
1.36025e-10 against unchanged tolerance 1e-10. The exact source, six successful
records, one rejected record and logs remain in `outputs/moving_current_initial_attempt`.
The final audit verifies its source archive; this attempt remains incomplete.

The final solver applies explicit linear iterative refinement to the original
reconstructed residual. One refinement is needed at nr128 and in the tight
nr64 control. Stable physical energy replay differentiates the same functional
before quadrature. Neither fix relaxes tolerances, changes physical equations,
alters density nor adds heat. All nineteen final cases use one frozen source.

## Reproduction and next gate

From the repository root:

```sh
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -m examples.moving_current_measurements
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -m examples.verify_moving_current
.venv/bin/python -m pytest -q
```

The driver accepts `--group moving_128` or another exact protocol name. Inputs,
kinematics/n_t, B_t/J_t/E, omission fields, configurations, versions, record/NPZ
hashes and source are archived. `outputs/moving_current_verification/report.json`
includes all audits; `checkpoint.json` records the inspected figure, tests and
source snapshot. Engine/reference/measurement fingerprint:

`d3984fb056edb1aa150d030937b9bce6c21d292b84ef80305b19b52f88af5489`.

Initial incomplete source:
`800d60a6a63899889e7c877b3f4c5bc10993cbb042b63b6f3b797a0260eb4012`.
Previous current/thermal source remains unchanged at
`80fd45bddda51871f7f36682e8d0bdc4495d0c41430191c25ef03e06ed47c4d2`;
historical records were not relabeled as new measurements.

Next implement **discrete conservative moving coupling**: actual continuity
flux supplies n_t; total mass is (mi+me)n; drift stress enters momentum and
its cylindrical terms; bulk/current/enthalpy/Poynting fluxes enter one total-
energy account. Verify local pressure partition, finite-mass compressible waves
and moving manufactured trajectories with time refinement and atomic rejection.
Do not insert this rate into the old split by changing only density.
All jobs are terminal. Stage2 and the full FRC goal remain open, including
transport/ionization, heated free edge/exact vacuum, real ends/coils, shock/
entropy and kinetic validity, formation and interactive experiments.
