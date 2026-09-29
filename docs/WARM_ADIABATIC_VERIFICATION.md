# Moving adiabatic species: measured verification

Status: **verified progress**, not a warm Hall or formation completion claim.
The full suite passes **283 tests in 94.53 s**. All **36 evolution records**
(30 main, 6 supplemental) complete and pass source/configuration/archive/history,
independent physical-energy and species-entropy replay. Both exported figures
were visually inspected. No numerical jobs remain active at this checkpoint.

## What was tested

The conservative mass overlap commutes **locally** with flux divergence,
including axis/wall half supports and the periodic seam. Independent affine
flux integrals and arbitrary open-boundary fluxes verify the mapping. Actual
warm RK stages reproduce the prior ideal fluid state bit for bit at both
orders; spatially constant specific entropy follows changing density and the
analytic gamma pe ndot/n compression response. Invalid final species states
leave all accepted fields and accounting untouched. Resistivity is rejected
until a physical species heating partition exists.

Independent references use fully filled R=.035 m, L=.24 m cylinders with
argon mass, n0=1e20 m^-3, Te0=5 eV, Ti0=1 eV, Bz=.01 T. These are verification
parameters, not validation of the experimental transport regime. Contact
translation has 25% density / 12% pe variation, v=2400 m/s, and evolves one
period. The finite-amplitude simple wave has velocity amplitude .05 cs and
evolves half an acoustic period before characteristic crossing. Reference
continuity, momentum, total energy and both species pressure equations pass
independent differential checks at gamma 5/3 and 1.4. Initialization integrates
conservative quantities with Gauss order 12, independently replayed at 18.

## Physical errors and representation

Table entries are continuous-volume relative L2 errors of the **pressure
variation**, not errors divided by the larger background pressure. Electron pe
is Q1 on the nodal grid. Reconstructed pi is reconstructed total pressure minus
that same pe; stored ion cell-mean errors are separately retained. All radial
cells are included. Contact electron/ion errors coincide because total pressure
is constant and species perturbations are opposite.

| nz | Contact pe / pi | Compression pe | Compression pi |
|---:|---:|---:|---:|
| 16 | 0.465559 | 0.108536 | 0.324027 |
| 32 | 0.132737 | 0.0270137 | 0.0900635 |
| 64 | 0.0449949 | 0.00717180 | 0.0229054 |
| 128 | 0.0169256 | 0.00201060 | 0.00586848 |
| 256 | **0.00561259** | **0.000592550** | **0.00150679** |

Finest Te/Ti perturbation errors are 0.00230161/0.00605575 for translation and
0.000670223/0.00343809 for compression. Finest node-only pe errors are
0.00567502 and 0.000579296; these are not substituted for continuous-field error.
Pa/eV RMS errors and positive reconstructed minima are in the report. Order8/12
physical-norm changes are at most **2.41243e-10** absolute in relative-L2 units.
Observed contact convergence is below uniform second order across this sequence;
the limiter and two-grid transport errors are measured, not hidden by the guide.

At nz128 first-order donor/Rusanov controls have node errors 0.891497 contact
and 0.0793516 compression, versus second-order 0.0172499/0.00193797. This is
declared numerical diffusion, not physical resistivity or collisional heat.

![Physical convergence](../outputs/warm_adiabatic_verification/warm_convergence.png)

## Time, stronger gradients and entropy defect

Fixed-grid nz32 controls use 128/256/512/2048 steps through .05 periods.
Electron-pressure differences from the 2048-step trajectory have ratios
4.03370/4.50765 contact and4.10340/4.30592 compression. Their coarsest RMS
differences are 2.12748e-5 Pa and 2.73992e-6 Pa. Finite reference-time error and
limiter switches affect ratios; these are separate from continuum spatial error.
At nz128, halving CFL changes pe by 0.000469704 Pa contact and 0.000383365 Pa
compression: 0.401%/4.19% of their respective base node errors.

The stronger contact has 60% density/16% pe variation. It exposes a negative
physical electron entropy change. Here the entropy column is
**(gamma-1) delta Se/(N kB)**, from the discrete logarithmic entropy, not the
conserved linear proxy U. No entropy inequality is claimed.

| nz | Strong contact physical pe error | Entropy change per particle |
|---:|---:|---:|
| 64 | 0.152234 | -5.83338e-4 |
| 128 | 0.0417378 | -1.09066e-4 |
| 256 | 0.0123756 | -2.39562e-5 |

The defect decreases but remains finite. Strong Mach .12 compression physical
pe error falls 0.0100661→0.00265388→0.000738578 at nz 64/128/256. Every accepted
species cell remains positive; the minimum reconstructed ion pressure across
all 36 cases is 2.22255 Pa. Stronger-gradient, near-vacuum or shock positivity and
species entropy robustness are still unproven.

Reference nref 1e20→7e18 changes pe by 1.06722e-15 relative. The nr 4→8 control
initially differs 5.14219e-6 because changing radial CFL changes the time step.
An additional **identical-time-step** control reduces the difference to
1.61383e-15. Opposite translation, gamma 1.4 and nearly cold electron
Te0=1e-6 eV controls also finish without floors or changes to fluid dynamics.

## Conservation and provenance

Maximum accepted-history mass/total-energy/proxy-integral relative drifts are
5.40031e-16 /8.87805e-16 /1.12490e-15. These roundoff budgets do not remove the
measured species partition or entropy errors. Independent direct magnetic
Gauss integration and explicit quarter-cell density/corner-pressure sums replay
thermal energy, both species energies and logarithmic entropy from raw fields.

Engine/reference/measurement fingerprint:
`0e37c39afe2dc22bd122cd6d21ec5e21981cfa97b7d5d59e154bd263d88c73b4`.
Main and supplemental JSON/NPZ collections are in
`outputs/warm_adiabatic_verification` and `outputs/warm_adiabatic_controls`.
The supplement driver and all report/audit/figure files have separate hashes.
Accepted histories use lossless scalar/nested columns. All 36 cases were
regenerated after compacting this layout; every accepted value, endpoint,
metric and NPZ hash matched the earlier layout exactly. Both inspected figure
hashes are unchanged. This reduced the saved collections from about 143 MB to
43 MB without dropping diagnostics; `compact_history_replay.json` records
that one-time comparison. The fingerprint above is the final regenerated one.
Reproduce with:

```sh
.venv/bin/python -m examples.warm_measurements
.venv/bin/python -m examples.warm_controls
.venv/bin/python -m examples.verify_warm_adiabatic
.venv/bin/python -m pytest -q
```

The prior 63 Hall records retain their original fingerprint. Since exposing
mass flux changes the toroidal module's source identity, their source hashes
were independently checked against commit 16969a9 rather than falsely relabeled
as current-source measurements. The isolated electron 46-case fingerprint is
unchanged. Existing Hall and electron regressions pass in the full suite.

## Remaining physical integration

This model uses common velocity and ideal MHD electric boundary data. Even
these exact one-dimensional references have a curl-free physical pressure
electric gradient that would be nonzero along a full conducting radial wall;
that full two-fluid boundary is **not** validated by this ideal fixture.
The pressure/current and Hall components therefore remain separate pending
consistent magnetic trace conditions and local conservative enthalpy/Poynting
fluxes. Warm Hall waves and nonlinear response, physical collision/heating
partition, conduction, ionization and full finite-mass transport still follow.

The original heated edge, true vacuum, physical gas ends, spatial coils, robust
shock/remapping/entropy behavior, regime validity, converged FRC formation and
dynamic experiment UI remain required. Stage2's formation guard stays in place.
