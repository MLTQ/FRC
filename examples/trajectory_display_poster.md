# Standalone scientific 3D poster

Reads the hash-verified float64 display arrays and full-vector magnetic paths
for step 17. Produces a PNG directly with Matplotlib, without requiring a
browser or advancing any physics. It is a complementary static view of the
same state shown in the interactive WebGL viewer.

Current density is shown on meridional cutaway faces, an outer cylindrical
surface and end faces, using actual |J| and a linear kA/m² color scale.
Surface colors are unshaded so they retain their quantitative legend. Surface
subsampling limits polygon count; it is only presentation. Magnetic lines use
all saved path points and are colored by total |B|, including the bias.

Physical aspect ratio is preserved, with one axial coordinate and two transverse
coordinates. The wire cylinder is a domain guide. The panel caption distinguishes
surface sampling from the interactive volume rendering. PNG annotations retain
the failed full-duration status and axisymmetric verification scope.

The output is `outputs/trajectory_display/current_state_3d.png`, separate from
browser-generated poster.png and all historical posters. A manifest captures
all inputs and output hashes. Run end-to-end, inspect the PNG, then freeze this
source. No unit tests are added or run.
