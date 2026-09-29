# Current-response case contracts

Defines 35 cases: uniform decays for four de values and
three meshes, one nr=64 check, twelve collisionless screening cases, three
variable-coefficient decays and six separate time/matching/angular/quadrature/
grading/mode controls, plus a finer nr=64 thin-screening case added after the
first measurements exposed spurious negative field on the unresolved nr=8 mesh.
The final complete protocol is rerun under a single frozen source fingerprint.
Complete typed inputs and physical interpretation are
recorded. Zero de explicitly means a massless model comparison.

All geometry is fixed. Decay initializes by nodal interpolation of analytic a,
so weak initial current error is measured, not silently removed. Screening
starts from zero and ramps a uniform analytic external field linearly. Variable
profiles follow the exact manufactured continuum reference, not a fitted
transport law. MaterialCurrent and all three analytic references share SI units.
