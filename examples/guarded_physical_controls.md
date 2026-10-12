# Explicit amplitude and residual controls

`run(destination, control, nr, nz, epsilon, tolerance, original=None)` executes
the complete verified component-stress signed workflow at F64/256 and fixed
radial16/N128 reference. Field counts, positive amplitude and original/tighter
solver tolerance are explicit, stored identically in physical and control
records. Tolerance cannot exceed the preceding 1e-10 bound. Radial field cells
must subdivide the fixed64 fluid elements. No equation or numerical operator
changes; the calculation retains its already verified physical source identity.

The separate control archive explicitly follows every dynamically selected
driver and physical source. It keeps 8 GiB own-process peak RSS and 256 MiB
free-disk guards, with SIGINT and five-second own-worker exit130 fallback.
Recorded final extrema must also satisfy both limits. Run one large solver or
independent audit at a time. Runtime/resource observations need not replay.

`--original` requires complete same-source protocol replay via the existing
all-array verifier. Accepted full signed rates, terminal strict rejections and
resource interruptions remain distinct. A complete rejection replay records
`accepted_physical_replay=false`; it does not create an accepted physical rate.
Incomplete resource interruptions cannot pass protocol replay. Every partial
attempt and exact input remains; fresh output paths are required.

The underlying driver saves inputs before solving and rates before q16/20
volume/face postprocessing. No result is rescued by a looser bound, filtering,
floors, modified masses, atmosphere or energy repair. Amplitude is a verification
perturbation, not a selected laboratory waveform or formation experiment.
Independent component/actual-curl comparisons and full replays are required.
No unit tests are added or run; old calculation sources remain frozen.
