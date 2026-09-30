# Preserved face-reconstruction failure

The first shared bulk-entropy implementation passed seven focused tests and all
14 instantaneous rate cases, but failed a full translating transit. This is a
numerical instability, not a failure of the physical model or a reason to add
density floors. Its source, complete outputs, exception and partial audit remain
under `outputs/bulk_thermal_*_initial_attempt`.

Frozen source:
`a90606b784ad3ec01fc3eda3897b60330bac775c8f217fbe6bf650e69100fbc5`.
The complete suite at that source passed 368 tests in 124.67 s. This did not
establish finite-time accuracy or stability of the new face interpolation.

## What failed

The initial normal-momentum faces reused the arbitrary-point density sampler.
That sampler averages two cell-containing polynomial stencils to preserve cell
means. At a normal face these stencils are not a symmetric pair around the face.
Floating-point floor decisions at nominal face coordinates can further change
which cell supplies the stencils. The resulting periodic translation operator
has a nonzero symmetric part and positive growth for one propagation direction.

The nr=8 full-transit density error was 0.249164 relative to the full physical
density. At nr=16 a later RK4 stage rejected nonpositive reconstructed density,
after the 0.75-transit snapshot. No floor or retry hid the rejection. The original
driver retained the traceback but not the last accepted state of a failed run;
that evidence limitation is explicit. The corrected driver now retains it.

Three compression cases completed, along with the inaccurate coarse translation.
Their initialization and last accepted step were replayed bit for bit, and all
saved snapshots independently integrated before the source was changed.
The nr=32 translation was deliberately interrupted after the instability was
independently diagnosed; the worker returned exit 130. Four remaining controls
had not started. The partial audit distinguishes all those states. It does not
present the initial ten-case protocol as complete or fully audited.

## Independent diagnosis and correction

`outputs/particle_face_diagnosis` contains 24 linear density-advection controls:
degree 3/5, three resolutions, both signs and both original/centered faces.
The complete operators, eigenvalues and full-transit density states are retained.
No thermal solver participates in this diagnosis.

For the original degree-5 faces, maximum positive growth exponents over one
transit are 12.94, 21.63 and 47.56 at nz=16/32/64. In the opposite direction the
real parts are near zero or negative. That directional asymmetry is retained;
the original negative-flow case is not mislabeled as unstable. Linear positive-
flow continuation at nz=32 produces negative density even while conserving total
particles. Global conservation cannot detect this error.

The correction uses one centered normal-face moment fit and integer face indices.
The periodic derivative is skew-symmetric for either constant-flow direction.
Its maximum real growth over a transit is below 4e-14 in the degree-5 controls,
consistent with floating-point error. No artificial diffusion was added. The
independent density-wave errors after a full transit are 1.61717e-4, 2.67701e-6
and 4.78902e-8 at nz=16/32/64, nearly identical for both signs. These are normalized
by the density perturbation, unlike the initial 0.249164 full-density error.

Independent tests also cover all-mode RK4 amplification at the declared step,
an arbitrary signal over forty transits using a matrix exponential, analytic
cell-average plane-wave convergence and radial regular moment fits. They do not
claim general nonlinear positivity, shock stability or stability of arbitrary
radial flow. The corrected full compression/translation protocol must still
pass its own physical and source audits.

The mean-preserving density reconstruction, physical entropy metric and
divergence-commuting flux primitive remain unchanged. Only the normal momentum
face trace changes. Initial pressure-force probes are preserved separately;
they are diagnostics, not an adopted full mechanical update.
