# hops_drive.py

## Purpose
HOPS means **Harmonically Oscillating Plasma Sphere** (user, 2026-10-22).
Causal electrode feedback only; it never sets particle forces, density, shape,
temperature or a harmonic well. Poisson determines the actual field.

## Contracts
`CageDrive.initial` gives the initial voltage. `advance` observes the old ion
state and advances only the electrode voltage to the next step. Legacy IEC
waveforms remain identical. New DC/fixed/HOPS modes start at zero voltage.
All modes use the same positive ceiling and first-order voltage response;
the actual voltage enters the electrode-work ledger.

HOPS starts on, observes a warmup, then freezes the measured reference RMS
radius. It switches off on inward crossing of the release fraction and on at
outward crossing of the return fraction. A low-pass radius sensor and minimum
dwell are explicit controller inputs, not smoothing of charge/field/particles.
The measured radial speed is about the cloud's measured centroid. No radius
power law or resonant frequency is imposed. Absorption can bias live-cloud
moments; remaining counts and loss channels must be reported. If the cloud
does not contract or re-expand, the controller does not invent a pulse cycle.

`fixed` uses a declared period/duty after the same warmup, not a measured
resonance. `dc` stays on. Observations act one step later, with no future data.
Saved switch times define actual intervals; at least two post-warmup onsets
are needed before reporting a measured repetition frequency.

The single-pulse response probe uses `fixed`, zero warmup and a period beyond
the physical endpoint. Release is evaluated at a step boundary; floating phase
arithmetic can delay it by one timestep. Report actual recorded command times.
The Oct 25 5 ns command/2 ns response fixture releases at 5.010/5.005 ns for
10/5 ps steps. Its voltage agrees with the independent exponential response
to those actual commands, rather than an assumed unquantized release time.
