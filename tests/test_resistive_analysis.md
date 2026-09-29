# Independent accuracy-comparison sampling checks

A globally linear a field and independent linear current/r profile have exact
values after an affine deformation of nested spherical sectors. The test
computes reference points from shell radii and trigonometric vertex formulas,
then checks the sampler's physical gradient, J and composite volume quadrature.
An independent revolved-boundary frustum integral supplies the volume target.

Separate explicit cell-index sums verify pressure volume restriction and
integrated physical heat restriction on nonuniform cell volumes. These checks
validate the comparison machinery, not the accuracy of a plasma simulation.
