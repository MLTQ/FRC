# Driven-field convergence analysis

Loads only completed pulse records whose saved-field hashes match. `compare`
uses exact cylindrical volume restriction for primitive cell fields on nested
meshes; `restrict_cells` checks preservation of their volume integrals. Equal
grids use the same operation with factor one. Reported L2 differences are
normalized by the finer/reference velocity, current or perturbation from the
declared initial rho=1, p=1 and uniform bias. They are inter-grid differences,
not errors relative to a known exact coil-driven MHD solution.

Nodal current is compared on coincident nodes using cylindrical dual weights,
with half weights at the material boundary. The axis has zero volume weight.
The pointwise current at the material edge has its own lower-order limitations,
so its measured norm must not be described as a cell-volume conservative field.

`summarize` retains peak-speed position, physical kinetic energy and density
range at each saved time. RMS current depth uses J² times cylindrical nodal
volume and the distance to the nearest gas boundary; it is a measured current
width, not a universal skin-depth formula. Numerical/physical eta is local and
separate from integrated numerical/physical heat. Inputs and complete fields
remain in the originating records rather than being replaced by report tables.

`compare_ends` uses coincident, equal-resolution cells in a fixed central axial
window when the actual gas half-length changes. It rejects incompatible grids
and records the physical window. Whole-domain source work is compared separately;
neither that quantity nor the window implies an unchanged physical experiment.
