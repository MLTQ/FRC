# Full finite-current kinetic slab evolution

## Purpose

Run the original20ns,±10mm,1e14m^-3 physical Ar+/electron slab with Ti/Te=.2/2eV,
zero initial current and field, remote ends±.2m. The collisionless kinetic model
evolves current, parallel heat and pressure anisotropy. This is independent of
the fluid gamma5/3 closure; neither is an exact baseline for the other.

## Interface and outputs

`python -m examples.vlasov_slab FRESH_OUTPUT --cells 1280 --velocity-cells 64`
uses self-consistent Vlasov-Poisson. `--free` disables electric acceleration for
the independent exact characteristic experiment. `--sigma` sets the explicit
velocity interval (default8 thermal speeds); `--cfl` defaults.65. Species grid
uses its own physical thermal speed. No masses, duration, densities or data are
fitted to the comparison. Fresh destinations keep all failed attempts.

- `initial/final.npz`: full phase densities, grids, physical profiles, actual E.
- `history.npz`: every accepted dt, N/P/E, exact-P0 field energy, spatial and
  velocity outflow ledgers, and actual midpoint electric boundary-stress impulse.
- `profiles.npz`: all moments and E at1ns intervals, including t0 and full20ns.
- `exact_characteristics.npz`: independent analytic free-streaming cell moments.
- `summary.json`, source/runtime archive: inputs, completion/failure, measured
  conservation and free-transport errors. Source changes invalidate success.

The Gaussian tail probability and second moment omitted by finite velocity
limits are reported. Initial P0 velocity representation has a resolution-dependent
energy error relative to the exact truncated Gaussian, separately from temporal
energy drift. Nothing is renormalized. Final endpoint raw-moment relative L2
errors use full spatial support and all velocity moments0..3. Number comparison
uses the independent characteristic integral over the whole spatial domain.

Zero inflow and actual outflow implement absorbing computational boundaries,
not physical plasma-wall sheaths. Electric stress impulse is included in total
momentum balance. Energy balance includes lost particle energy; if charge exits
in a future case, external field boundary work must be derived before treating
that particle-corrected ledger as closed. Velocity-domain charge loss can also
change potential energy; its observed amount must be reported. Current cases
have remote ends and velocity losses below2e-16 of initial ion number.

## Resource and scientific contracts

At least384MiB free is required before launch. Every250 accepted steps checks a
320MiB export reserve; failure saves accepted state and never labels it20ns.
Do not start concurrent copies of expensive fine runs. No unit tests are run.
This driver assesses a new physical component, not a verification framework.
Explicit phase-space/time and velocity-domain sensitivity remain necessary.

| Depends on | Expects |
| --- | --- |
| `plasma_bench/vlasov_slab.py` | conservative positive transport; no state repair |
| `plasma_bench/ballistic_slab.py` | independently integrated raw cell moments |
| `implicit_plasma_records.py` | source/runtime identity and preserved sources |
