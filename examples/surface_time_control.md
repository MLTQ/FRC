# Supplemental finest-mesh time control

Reproduces the additional nr=32, grading=2 run with CFL=.0375 and explicit
time cap .000125 s. It was first launched by a direct call to the existing
frozen `measure(ResistiveCase(...), ..., fingerprint=provenance)` entry point,
after the initial fine control showed noticeable time sensitivity. The actual
typed configuration is completely recorded and checked by the assembler.

The production/shared-measurement fingerprint is the same as the original
20-case study; this later reproduction helper is separately hashed with the
report. No physics or measurement algorithm changed to add this case. Run
`python -m examples.surface_time_control` only when intentionally reproducing
it; do not launch a duplicate while its original job is live.
