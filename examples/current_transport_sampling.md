# Independent physical transport reconstruction

Reconstructs compatible current and magnetic fields by root-product polynomials,
including second derivatives of nodal bases for edge-current gradients. Regular
fluid scalars and dependent transverse axis traces use separate root polynomials;
periodic axial fields use complex FFT coefficients with the real Nyquist cosine.
No production sample, gradient, adjoint or mass solve is called.

`TransportSampler` provides physical cylindrical quadrature, fields, covariant
vector gradients, and the exact derivative of the declared velocity quotient
for a supplied state rate. It checks positive sampled density and real periodic
reconstruction. Current field degree is fixed at five, matching this study's
protocol; fluid degree is read from the fixture. Field cells subdivide fluid
elements, as required by the original mapping. Integration excludes r=0 without
a radius floor. Independent analytic targets provide exact axis limits.

Periodic points are wrapped before stencil selection; translations by a full
period therefore retain the same polynomial branch. A focused test exposed the
initial out-of-domain upper-cell clamp; that audit-only defect was corrected
before freezing this study. No production transport source changed.
