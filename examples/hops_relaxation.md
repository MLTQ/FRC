# hops_relaxation.py

## Purpose
Read-only measurement of a single-pulse HOPS ion response from every archived
native particle state. Distinguishes live-ion survival in the chamber from
occupancy inside the open wire cage. Does not evolve or alter the simulation.

## Components
- `measure`: uses full native ion positions and velocities; measures cage
  occupancy, core counts, radii, outward motion and actual potential history.
  Checks ion absorption-count closure and recomputed radius/speed against history.
- `first_crossing`: reports first sampled threshold and preceding timestamp,
  without interpolating a state. Thresholds not reached have a finite observation
  window, not an invented exponential fit or infinite lifetime. Conclusions use
  measured dilation signs and threshold results, not an assumed expanding cloud.
- `compare_archives`: optional bitwise comparison of every saved array at all
  common original archive timestamps. This measures archival independence only.
- `plot`: creates standalone PNG/PDF scientific figures from the measured CSV.

## Definitions and contracts
- Cage occupancy uses the convex hull of the wire centreline vertices, whose
  faces are open in the physical solver. This diagnostic adds no reflecting or
  absorbing surface. A separate sphere at centreline circumradius + wire radius
  gives an unambiguous count beyond the entire cage.
- Ion RMS radius is centroid-relative. Radial quantiles and core counts use
  the chamber origin. Absorption biases live-population moments; these are not
  individual particle tracks or a fixed survivor cohort.
- Rdot=<dx·dv>/R measures instantaneous kinematic dilation, not the derivative
  of a lossy population across absorption events. No symmetry is imposed.
- R/Rdot is an instantaneous size-change scale, not a damping/e-folding time.
  The fitted exponential is solely the prescribed electrode voltage fall.
- Inputs need history, geometry, drive switches and all native state files.
  Assumes equal macro weights, no ion source, a single off edge and nonempty ions,
  as in the declared pulse-release experiment. Output directory must be new.
- Outputs: `time_series.csv`, `relaxation.json`, `relaxation.png/.pdf` and optional
  `shared_times.json`. Source files and scientific inputs are preserved.

## Use and limits
`.venv/bin/python -m examples.hops_relaxation
outputs/lab_runs/hops_pulse_release_1ns_v1 outputs/hops_relaxation_20261101_v1
--reference outputs/lab_runs/hops_pulse_release_v1`

No unit tests or solver rerun are needed for this analysis. No collision or
thermalization time can be measured from a collisionless solver. Coarse-grained
collisionless damping could still occur, but it needs an observed response and
appropriate analysis; no bounded post-release cycle is assumed. Spatial-force,
particle, phase and wire-interception accuracy remain unresolved.
