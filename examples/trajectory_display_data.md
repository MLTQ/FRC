# Saved trajectory field display data

## Purpose
Read the fully replayed refined contact trajectory and reconstruct four actual
states for a 3D display. This read-only postprocessor never advances or changes
the simulation. It preserves the failed full-duration status.

## Components
- `FieldSampler`: independent degree-five root-polynomial reconstruction of all
  magnetic and weak-Ampere current components on arbitrary tensor grids.
- `trace`: cylindrical full-vector field lines, using z as parameter only after
  checking positive axial B. RK4 follows the interpolated display field; the
  half-step endpoint comparison measures tracing error, not physical accuracy.
- `run`: validates the complete original/repeated trajectory, computes current
  with the unchanged physical operator, compares independent samples against
  the production quadrature, archives float64 display fields and compact JSON.

## Contracts
The source is the F8/64, B32/64 contact fixture, radius 35 mm, length 240 mm,
periodic z, filled plasma, 5 mT bias. The last accepted state is step 17;
step 18 rejected. There is no FRC formation, vacuum edge or 3D stability claim.
Four frames are actual saved steps 0, 10, 15, 17. Common linear color scales
cover all frames. Current opacity and line seeds are presentation choices.
No field-line deformation, filtering, numerical diffusion or state repair.

The grid uses radial and axial cell centers. Samples outside that display
support are excluded from interpolation audits and traces; rendering clamps
the thin outer half-pixel strip explicitly as a presentation approximation.
Uniform physical length scaling preserves the cylinder aspect ratio.
Both signs of every vector component are retained in the raw NPZ files, while
volume color uses the Euclidean magnitude in A/m². B includes the axial bias.

Independent sampling versus production and off-grid interpolation errors are
reported separately. Off-grid comparisons use actual field differences on a
fixed interleaved tensor grid; they do not certify physics convergence.
The complete postprocessor is replayed into a fresh directory and compared.
Once executed, preserve this source; subsequent scientific variants get new names.
