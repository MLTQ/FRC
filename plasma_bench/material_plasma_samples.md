# Interior fields at common material labels

`material_samples` evaluates the actual nodal/histopolant magnetic and electric
bases of MaterialVacuumSpace at arbitrary strictly interior reference r0,z0.
The actual Q1 coordinate/angle maps give F, B=F*B0/det(F), and w=F^-T*w0.
The weak-Ampere drift coefficients are reconstructed from the saved state, not
replaced with a finite-difference curl of a displayed field. J=e*n*w with
n=n0/det(F); the caller must supply its own declared initial density function.

Flux is the actual reference poloidal polynomial plus the exact prescribed
source flux at the moved point. Physical velocity uses actual nodal velocity
and angular speed. Coil fields use the source waveform at the requested time.
No field evolution, cloud interpolation, smoothing or state correction occurs.
This supports actual field plots and material-label spatial comparisons between
different meshes, including complete B/J/flow components and moving coordinates.

Inputs exclude the axis and exact plasma boundary, where separate regular/trace
limits are necessary. Element interfaces select the right cell; finite-space
derivative jumps remain. Identical material labels on two runs can map to different
physical positions: report those differences alongside field/current norms.
Use an independent common quadrature, explicitly weighting its reference volume
by a chosen run's actual determinant. Such differences are not experimental
accuracy or proof of convergence from only two grids.

Requires MaterialVacuumSpace coefficient ordering and the complete saved fluid
angle. The compact MaterialVectorSpace has a different poloidal representation
and is not supported. Callers may reuse an actual `response` to avoid another
weak-Ampere reconstruction. These field samples are not particle trajectories
or topology certificates.
