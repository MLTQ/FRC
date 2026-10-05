# Exact one-sided physical traces

Independently builds root-product field interpolation at the exact supplied
coordinate, with an explicitly selected field cell for each point. Fluid
radial-element selection follows aligned mesh topology, not a displaced floating
point or an average. Left/right traces can therefore differ without ambiguity.
Periodic axial cells may use index -1 at the first seam, with translated roots
and wrapped coefficient indices. Fluid axial sampling retains all Fourier modes.

Inherits only physical reconstruction formulas from TransportSampler; its
sampling matrices are built independently for exact trace locations. Supports
actual density/velocity/entropy, current/electric and magnetic coefficients.
No state or flux is modified. The axis is excluded because its cylindrical face
area is exactly zero; the calling face audit inserts that zero geometrically.
The wall uses its interior cell. Matching radial fluid-element/field-cell
alignment is required. Manufactured smooth and discontinuous trace tests are
independent of any evolving plasma or fitted mode.

Fluid cardinal polynomials use direct products at the requested point rather
than evaluating expanded powers near an endpoint. This preserves constant
trace values to roundoff without imposing a partition-of-unity correction.
