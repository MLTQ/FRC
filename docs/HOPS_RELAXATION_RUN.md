# HOPS single-pulse release: 25 µs experiment design

Status: **designed and waveform-validated; production cases not launched**.
Tracked by beads `frc-cac.2.7`. This tests the slow ion expansion seen in the
completed 5 µs pulse experiment before choosing adaptive HOPS timing.

## Question and hypothesis

Does the released cloud keep emptying through the open cage, or does its
self-consistent field produce a coherent return/compression phase? Measure the
response directly; do not assume an exponential relaxation, harmonic well or
radius-frequency scaling. A long-lived population alone is not an oscillation.

## Two physical cases

| Input | Pulse case (first priority) | Matched control |
|---|---|---|
| Actual physical interval | 0–25 µs | 0–25 µs |
| Cage command | +3 kV from 0 to 0.2 µs, then zero | Zero throughout |
| Voltage response | 20 ns first-order rise/fall | Same operator, zero target |
| Ideal wall electron emitter | 1 mA, 2 eV, stops at 0.2 µs | Identical |
| Initial plasma | Same seed 41, 600 macros/species | Identical |
| Initial density/radius | 1e11 m^-3, 5 cm | Identical |
| Species | Actual-mass e and singly charged mass-2mp ions | Identical |
| Initial Te / Ti | 1 / 0.05 eV | Identical |
| Cage/chamber | 15 cm open wire cage / 40 cm grounded sphere | Identical |
| Mesh | Existing 18 mm bulk / 1 mm wire mesh | Identical |
| Integration / native archives | 5 ps / 1 ns | Identical |

The control removes the electrode pulse while retaining the initial plasma and
prescribed source. It measures the pulse's effect within this idealized model;
the two cases do not have equal delivered work. Compare their actual ledgers.
It is not an emitter-free control or the historical air discharge.

Use existing fixed-drive support with **40 µs period and 0.5% duty**. The product
is a 200 ns pulse; the next pulse is outside the entire 25 µs interval. The old
10 µs period must not be reused for this experiment. The 40 µs setting is only
a way to issue one pulse, not a proposed HOPS repetition period. Source cutoff
is integration step 40,000; there is no injection or second pulse afterward.

## Why 25 µs

The 5 µs observation ended before cage half-empty or e-fold occupancy. Starting
from its native 0.5 µs ion state, an independent no-force straight-flight
reference reaches half the initial cage occupancy near 16.57 µs and 1/e near
20.25 µs. This motivates a 25 µs observation window. It is a duration guide,
not a prediction that the self-consistent plasma must follow those times.

For that reference, dx and dv are relative to the measured centroid/mean
velocity, and

`R²(t) = <dx²> + 2 <dx·dv> Δt + <dv²> Δt²`.

It deliberately omits both forces and metal absorption. Do not mistake it for
a physical loss control. Later live-population RMS values can differ simply
because the simulation absorbs its fastest particles.

## Measurements and interpretation declared before execution

- Cage occupancy: count inside the wire-centreline convex hull, normalized to
  the original 600 ions. Report first 90%, 50% and 1/e crossings and whether a
  crossing reverses later. Hull faces are diagnostic only, physically open.
- Retention: live ions anywhere in the chamber, live ions outside the cage,
  and cage/wall absorption separately. Also count beyond the full 7.5 cm outer
  cage radius. Show physical counts/densities using the saved macro weight.
- Compression: fixed-centre core number/density within 2 cm, centroid-relative
  RMS and principal radii, origin radius quantiles and centroid translation.
  A shrinking live RMS caused by selective absorption is not compression.
- Phase: ion radius versus radial velocity from all native ions, mean dilation
  Rdot and outward-moving fraction. Candidate return requires an inward phase
  supported by core-density and phase-space evidence, with loss changes shown.
  Inspect potential/E structure rather than inferring a well from one scalar.
- Timing: measure crossings from the actual 0.2 µs off edge. A period requires
  repeated equivalent turning points; fit damping only with a resolved sequence
  of cycles. If motion stays outward, report expansion/emptying times. If a
  threshold is not reached at 25 µs, report the finite 24.8 µs release window.
- Power/accuracy: retain raw source/electrode work, field and species kinetic
  energy, injection, absorption and charge balance. Also measure incremental
  energy-ledger drift relative to 0.5 µs, normalized by live K+space-charge
  energy then. More than 1% late drift makes quantitative late dynamics
  unresolved and calls for timestep refinement; never repair the state.

These are observations, not feedback constraints. Run to the declared endpoint
without releasing/re-pulsing at a chosen radius. Empty populations are a loss
endpoint, not a zero-radius oscillator.

## Execution and fidelity

Configs: [pulse](../configs/hops_release_25us_pulse.json) and
[zero-bias control](../configs/hops_release_25us_zero_bias.json). Both are valid
inputs to the existing `examples.iec_3d` CLI with the existing mesh
`outputs/lab_runs/hops_pulse_release_1ns_v1/mesh.npz`. Freeze source and configs
using the existing campaign workflow, then run cases sequentially with one
compute thread and detached logs. New identifiers should be
`hops_release_25us_pulse_v1` and `hops_release_25us_zero_bias_v1` under
`outputs/lab_runs/`; do not overwrite existing runs.

The checkpoint loader binds the exact configuration, including duration.
This design therefore uses a separate from-zero replay with the existing
driver. Do not edit checkpoint metadata to pretend it is an unchanged resume.
An explicitly verified checkpoint branch could save the repeated prefix later;
it is not implemented or assumed here.

The changed period/duty preserve the original command prefix through 5 µs and
the actual dt is bitwise identical. During execution compare native particle
and field arrays at shared original timestamps before interpreting the new
tail. Waveform-only validation is not a substitute for that full replay.

Each case has 5,000,000 steps and 25,001 native states. The completed dense
5 µs run took 9.23 h and retained 2.56 GB of frame/state files. Linear planning
estimates are about **46 h and 13 GB per case**, with substantial dependence
on source-control populations and archive growth. About 41 GB was free at
design time. Preserve complete checkpoints and original outputs; no other
experiment is paused by this design.

Validation and reference CSV: `outputs/hops_release_design_20261101_v1/`.
The actual controller's first 0.5 µs voltage matches the old pulse bitwise and
the independent exponential within 1.1e-11 V. Full-interval gate arithmetic has
exactly one off edge at 0.2 µs and no subsequent on edge. This checks the
experimental command, not physical convergence. No solver change or unit suite.

## Scope

Evolve full 3D collisionless Vlasov–Poisson: dx/dt=v, m dv/dt=qE,
E=-grad(phi), -epsilon0 Laplacian(phi)=rho, with actual absorbing steel and
prescribed electrode potentials. No imposed symmetry, well, return flow,
neutralizing atmosphere, mass scaling or energy correction. Tetrahedral P1
force approximation, particle noise/self-force and wire interception remain
unconverged; numerical conservation does not establish their accuracy. There
is no neutral/Coulomb collision operator, chemistry, ionization or fusion.
This run tests collisionless expansion/return, not collisional thermal relaxation.
Any apparent resonance must survive timestep/mesh/particle comparisons before
it becomes a physical claim or a basis for choosing adaptive drive frequency.
