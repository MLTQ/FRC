# iec_cage_check.py

## Purpose
Actual 90-wire cage mesh refinement at bulk/wire targets 24/1.5, 18/1 and
14/0.75 mm, compared with the earlier independent boundary-method result.
This measures a concrete new solver's static error, not a recursive framework.

## Contracts
`python -m examples.iec_cage_check NEW_OUT [--reference EXISTING_18MM_MESH]`
saves each full mesh and report.json with capacitance and centre potential.
The independent MFS values are from `outputs/iec_vacuum_reference_v4`; they have
up to 1.4% local boundary potential error, so agreement is cross-method evidence,
not an exact near-wire gradient or interception baseline. The analytic sphere
and capacitor comparisons are in `iec_field_check`. Actual kinetic results still
need particle/mesh/dt convergence before tuning or confinement claims.
