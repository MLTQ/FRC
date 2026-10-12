# Directional field and memory report

`run(destination, cases, resources)` renders standalone scientific PNG/PDF from
complete signed calculations, independent current/strong-curl audits and their
full replay records. `cases` entries supply label, calculation directory, audit
directory and simulation verification JSON. `resources` supplies labels and
primary/replay guard-summary paths. Fresh report paths are mandatory.

Verifies terminal states, complete bitwise replays, matching audit inputs and
current input hashes, original physical fingerprints, q20 fields, untriggered
resource guards and final peak/free-disk limits. Archives source, input hashes,
plotted values and rendered file hashes. Inspect PNG visually after rendering.

Field-grid labels state radial / axial cell counts; plotted errors use q20.
Current and magnetic rates use each component's own continuum norm. Hollow
current markers are independent inside-cell curl(Bt)/mu0 and exclude face
distributions. Memory covers entire signed workflows and audits/verifiers;
run-to-run observations are not guaranteed bounds or exclusive allocations.

The fixture keeps fluid64/256, radial16/N128 reference and original1e-10 bound.
This report illustrates separate field directions, not asymptotic convergence,
full time conservation, moving vacuum, formation or laboratory suitability.
No numerical calculation, physical field or solver acceptance is changed, and
no unit tests are added or run.
