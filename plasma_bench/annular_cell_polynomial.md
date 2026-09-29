# Mean-preserving transverse polynomials

Radial and axial cell polynomials average mirrored four-cell cubic fits. Each
fit contains the target cell and reproduces its actual volume average; physical
interior stencils are clipped at the radial boundaries without adding exterior
matter. Radial moments use r^s P3(r²/2), s=0/1/2, with exact physical quadrature.

These are distinct from shared face interpolation in annular_polynomial: when
integrated transversely across a face, the polynomial must recover that face's
input mean. Arbitrary-point face-stencil switching would violate that property.
The operators are linear and preserve inputs; no limiter or positivity repair
is included. Smooth bulk transport validates positive reconstructed states.
