# Measured linear Hall and electron-inertia accuracy

The separate bounded axisymmetric component now evolves poloidal/toroidal field
and flow, counts exact finite-species-mass current kinetic energy, and measures
positive resistive heat. **This is progress toward the active formation goal,
not a completed FRC experiment solver.** Governing equations, restricted boundary
traces, FE/Fourier discretization, omitted physics and independent references
are in [LINEAR_HALL_NUMERICS](LINEAR_HALL_NUMERICS.md).

The full suite passes **208 tests (78.87 s)**, including 16 new checks. New tests
cover original species momentum equations, both helicities and bias signs,
cylindrical curl/divergence and regular axis, massless and finite-mass limits,
physical energy/heat integrals, angular momentum, spatial/time convergence and
atomic caller-state preservation. A purely poloidal initial perturbation produces
toroidal magnetic field and flow through Hall coupling; the MHD control does not.

## Complete records

`outputs/hall_verification` contains 49 complete cylindrical cases, per-step
histories and lossless initial/final coefficients plus sampled physical fields.
Independent parallel verification contains 16 modes, each with 40/80/160 steps.
All complete configurations and relevant source/archive hashes are checked by
`examples.verify_hall`. No smoothing, resistivity floor, density floor, anomalous
resistance or energy correction is used. Frozen cylinder measurement SHA256:

`868ece2de1730a51d4aaf363652e3a19139c8289e046eda3188b57a0a45b0e3f`

The parallel and report/figure sources have separate hashes. Historical current
and moving-plasma studies are preserved; their findings are not overwritten.

```sh
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl XDG_CACHE_HOME=/private/tmp/plasma-bench-cache .venv/bin/python -m examples.hall_wave_measurements
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl XDG_CACHE_HOME=/private/tmp/plasma-bench-cache .venv/bin/python -m examples.hall_parallel_measurements
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl XDG_CACHE_HOME=/private/tmp/plasma-bench-cache .venv/bin/python -m examples.verify_hall
```

## Cylinder convergence

Uniform n=1e18 m^-3, physical proton/electron masses, B0=.001 T and prescribed
eta=3e-5 ohm m are verification inputs. MHD and Hall fixtures use R=1 m, L=2 m;
the finite-mass fixture uses R=.01 m, L=.02 m to reach electron scales.
These are explicitly distinct mathematical references, not a substitution for
the original low-pressure argon formation target. All start from interpolated
continuum Bessel waves and evolve through two radians of their exact frequency.
The finite-mass reference is incompressible; it does not claim cold compressible
oblique-wave accuracy or kinetic applicability.

| Model / curl helicity | nr=64 field relative L2 | Flow relative L2 | Weak-current relative L2 | Relative phase-frequency error |
|---|---:|---:|---:|---:|
| MHD, either | 0.00625922 | 0.00625922 | 0.00625945 | -1.13475e-6 |
| Hall, h=-1 | 0.00625973 | 0.00625916 | 0.00625948 | -2.74641e-5 |
| Hall, h=+1 | 0.00625919 | 0.00626011 | 0.00625990 | 1.90887e-5 |
| Finite mass, h=-1 | 0.00626052 | 0.00625973 | 0.00625962 | -4.01628e-5 |
| Finite mass, h=+1 | 0.00625948 | 0.00625931 | 0.00626125 | -3.05562e-5 |

Field errors decrease from approximately 5.2% at nr=8 to 0.626% at nr=64:
first-order derivative-field accuracy for these piecewise-linear potentials.
The smooth phase errors approximately follow second-order radial refinement
until time error matters. Higher radial mode 2 / axial mode 3 cases also converge:
Hall h=-1 field errors are .046402, .022357, .011070 at nr=16/32/64.
Doubling axial samples 9→17 preserves the resolved single-mode result; this is a
Fourier representation check, not evidence that arbitrary data avoid aliasing.
Radial operator quadrature 8→12 changes the result negligibly.

Phase and damping come from an evolved complex amplitude projection. Interpolation
excites additional spatial/helicity modes: even a collisionless run can show a
small apparent damping or growth of this projection while total energy is
unchanged. Such mode interference is not physical collisional heating or a
negative-resistivity claim. Full-field errors and heat are reported separately.
At nr=64 finite-mass h=-1/+1, apparent damping error divided by reference
frequency is 1.74340e-5 / 6.09553e-8. The true damping rates are 460.275 and
738711.745 s^-1 respectively. These normalizations must not be misread as
relative errors in a small damping rate.

## Time, dispersion and energy

Independent parallel circular-wave tests span k di=.1/1/10/1000, both helicities,
with and without electron inertia. The actual time-halving error ratios are
3.998876–3.999721. At k di=1000, the Hall fast branch has omega=9.57884e10 s^-1;
physical finite electron mass limits it to 1.75560e8 s^-1. The corresponding
ion branches are about 9.57860e4 s^-1. Direct species-equation tests and separate
large-k tests recover the physical electron/ion gyrofrequency asymptotes.

The slow cylinder branch initially exposes unresolved fast secondary modes:
128→256 steps changes apparent damping nonmonotonically. A separately assembled
helicity-space dense exponential retains all initial spatial modes and isolates
time error. With 1024/2048/4096/8192 steps, its relative complete-energy norm
errors are 6.35558e-7, 1.58885e-7, 3.97212e-8 and 9.93030e-9. Thus the finer
controls resolve second-order time convergence without mistaking the remaining
spatial/projection error for a time failure.

Maximum measured whole-run energy-plus-heat residual is 6.63520e-15 relative;
maximum equation backward error is 2.38842e-17. Independent physical-volume
integration agrees with the modal energy within 5.12417e-14 of initial energy;
raising verification quadrature changes field norms by at most 3.09475e-15.
These energy checks do not establish phase accuracy: an independent one-step
collisionless test advances two radians with a **21.46% phase error**, while
preserving energy to roundoff. Finer steps are required despite stability.

![Measured Hall accuracy](../outputs/hall_verification/hall_accuracy.png)

![Actual toroidal wave and rotation](../outputs/hall_verification/hall_wave.png)

## What remains

The compressible/moving solvers still lack these toroidal degrees of freedom and
the nonlinear Hall/current-energy coupling. Uniform-density linear projection
cannot simply be transplanted into a compressing plasma. Next extend the bounded
compressible state with toroidal field and conservative axial angular momentum,
verify cylindrical stress/induction work and nonlinear ideal transport, then
add Hall and finite-mass terms with their parent-species work pairing. Separate
ion/electron thermal energies and edge transport remain necessary for formation.

No new result resolves the heated sharp-edge incompatibility, actual spatial
coil force, physical vacuum/end coupling, robust shocks/remapping, kinetic regime
or converged formation. The Stage 2 formation guard remains in force, and the
overarching goal remains active. An axisymmetric wave says nothing about 3D
stability or particle confinement.
