# Resistive material specifications and reproducibility

`ResistiveCase` declares ideal-initialized FRC relaxation and uniform Gaussian
bias-ramp fixtures, with physical scalar D and eta=mu0 D, full material/magnetic
geometry, outer timestep and adaptive diffusion tolerance. These unit-scale
inputs have no inferred species, temperature or transport validity.
Optional `radial_grading` and `angular_cells` expose independent surface-normal
and tangential resolution controls; defaults retain the uniform 2*nr mesh.
Configuration records always include the actual angular count and grading.

`snapshot` extends complete material archives with weak current Jtheta/r,
current quadrature positions/values and cumulative per-cell Ohmic heat. Entropy
and flux labels are saved at every archived time because they now evolve.
`provenance` freezes all production modules and actual measurement/reference/
snapshot helpers with Python/NumPy/SciPy/platform versions. Historical ideal
records retain their own source fingerprint; they are not relabeled or rerun
silently after the resistive extension.
