# Conditional physical transport and carrier checks

The experiment's prescribed conductivity does not become a species closure
merely because an energy budget closes. The new diagnostics evaluate separate
local Coulomb estimates from explicit ne, Te, Ti, mi and integer Z. They do not
overwrite conductivity or establish ionization equilibrium.

## Classical coefficient convention

The primary reference is the [NRL Plasma Formulary, 2009](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pages 28–29, 34 and 37–39. Its cgs/eV expressions are converted to SI
using physical constants; independent tests retain its rounded cgs coefficients.
For slow thermal ions and Maxwellian electrons, the implemented electron-ion
Coulomb-log branches use electron density in cm^-3:

```
lnLambda_ei = 23 - ln(sqrt(ne_cm) Z Te_eV^(-3/2)), Te < 10 Z² eV,
           = 24 - ln(sqrt(ne_cm)/Te_eV),          Te >= 10 Z² eV.
```

The required ion-temperature ordering is checked. Same-species ion-ion log is
`23-ln[Z³ sqrt(2 ni_cm) Ti_eV^(-3/2)]`. Nonpositive logs and unsupported
branches are unavailable; no arbitrary minimum logarithm replaces them.
With C=e^4/(4 pi epsilon0)^2,

```
nu_ei = 4 sqrt(2 pi) ne Z C lnLambda_ei / [3 sqrt(me) (e Te_eV)^(3/2)],
nu_ii = 4 sqrt(pi) ni Z^4 C lnLambda_ii / [3 sqrt(mi) (e Ti_eV)^(3/2)],
eta_momentum = me nu_ei/(ne e²),    lambda_thermal = sqrt(kT/m)/nu.
```

For Z=1, the separately named parallel/unmagnetized classical resistivity is
eta_momentum/1.96. The momentum coefficient is the perpendicular friction
limit for strong magnetization, not a complete scalar conductivity or the
Pedersen response. The factor is not extended to other Z. Near weak coupling,
large drift or non-Maxwellian distributions these expressions are unreliable.
Neutral collisions, thermal forces, ionization, radiation and transport tensors
are still absent; partially ionized records say so explicitly.

## Current-carrying and local validity scales

Derived directly from plasma frequency and relative current velocity,

```
de = sqrt[me/(mu0 ne e²)],
di = sqrt[mi/(mu0 ni Z² e²)],
lambda_De = sqrt[epsilon0 kTe/(ne e²)],
v_drift = |J|/(ne e).
```

Compare these to actual structure/gradient length and time, not just the global
chamber radius. Electron drift is also compared to electron/ion thermal speeds
and c. Magnetic-null regions cannot be justified by a remote bias magnitude.
Mean-free-path and collision-time checks are conditional Coulomb estimates.
Anisotropic magnetized transport needs directional length scales; the single
supplied length is a screening diagnostic, not a complete transport tensor test.

Ratios >=0.1, Coulomb logarithms <5 and light-transit/time >=0.01 trigger review.
Those are explicit diagnostic thresholds, not sharp physical transitions or a
proof of validity when no flag appears. `validity_established` remains false.
Vacuum has no transport state; zero-density input is rejected instead of filled.

Stage 0/1 output now includes a conditional reference using input temperatures,
plasma radius, bias, shortest waveform scale and zero initial current. The
historical inferred-Drude scales and prescribed conductivity remain distinct.
Local evolving current and gradient scales require separate calls to `regime`.

## What this enables, and what remains

Explicit conditional scans distinguish collisional, Hall, electron-inertia and
carrier-drift requirements before selecting a dynamic model. They do not choose
experimental density, prove a pre-ionization state, or justify removing edge
heat. The separate current-inertia component verifies one physical response
term, with drift kinetic energy and a genuine exterior vacuum. Hall feedback,
electron/ion thermal evolution, physical gas ends and a justified edge/kinetic
treatment remain required before a formation experiment is credible.
