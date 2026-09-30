# Common-grid physical comparison

Independent root-product interpolation of n,mz,Ce,Ci,mr/r,L/r² in each element
gives regular physical fields at arbitrary positive radii. Fourier resampling
puts different axial grids on the same uniform physical points. No production
derivative or flux defines the comparison. Sampled density and pressure must
remain positive and finite; no clipping is allowed.

`grid` uses high-order cylindrical r dr dz quadrature on a declared reference
partition. `differences` compares both species and all velocities in common
physical norms, including a fixed physical central region. Velocity errors use
background sound speed, including when a reference component vanishes. Other
fields use the reference norm. A finer numerical solution is a convergence
reference, not an analytic truth; manufactured exact comparisons remain separate.
