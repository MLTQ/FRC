# Resource accounting before larger field calculations

Builds the exact original compatible/wall geometry and normalized lift on a
short32-cell axial domain with256/512 radial cells, then inventories retained
sparse arrays and actual short-domain electric-mass LU factors. Original
DirectMetricWallInduction lift expressions are reproduced only for accounting;
no physical current model, preconditioner or derivative is accepted.

Sparse storage is extrapolated linearly in axial cells beyond the fixed
stencil width. Full uniform-metric storage is computed from the actual source
formula nz*(2*nr-2)^2 complex128 entries. Scalar, vector and gradient quadrature
arrays are separately counted. Larger-grid sparse LU fill and live nonlinear
temporary arrays are explicitly excluded: this is not a peak-memory bound.

The reported macOS ru_maxrss is the cumulative own-process maximum in bytes.
No other process or global memory is inspected. Source/data remain archived;
no physical or numerical setting is changed. Complete guarded physical runs
and fresh replays must separately establish actual feasibility and accuracy.
