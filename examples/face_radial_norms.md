# Independent physical norms for face radial momentum

Uses the independent root-product primitive for interior radial faces and the
independent density/regular-angular maps for cell quantities. It does not use
production kinematics or primitive matrices. Physical velocity, both pressures,
temperatures and total energy follow their actual equations. The inherited
measurement records full-volume/first-annulus/max errors and initial projection.

`rates` independently differentiates true K and both thermal energies using
actual face/cell momentum rates and n_t. Radial, azimuthal and axial acceleration
errors are separate, preventing small whole-domain norms from hiding an axis
failure. This sampler can audit saved state/rate arrays at higher quadrature.
