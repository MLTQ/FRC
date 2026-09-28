# Numerical operators
`field_from_flux` provides second-order nodal display/probe fields with a regular
axis extrapolation. `cell_divergence` uses integrated face flux differences;
telescoping guarantees div B to roundoff, including cells adjacent to r=0.
`stiffness` builds a symmetric positive definite interior energy matrix with
homogeneous boundary flux. `expand` pads an interior vector with zero boundary
values. The stationary solver depends on this ordering and energy convention.
No artificial physics or smoothing is included.
