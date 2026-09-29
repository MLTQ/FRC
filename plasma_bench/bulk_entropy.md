# Bulk electron entropy transport

`entropy_rate` advects U=n w with exactly the mapped bulk particle flux Gamma:
Udot=-div(Gamma w_upwind). Constant w therefore follows the changing density
without an independent continuity solve. Gamma includes the primal solver's
Rusanov mass diffusion. This common-velocity operator has no Hall/current drift,
pressure battery, collision, conduction, resistivity or electron inertia.

Order 1 is donor-cell. Order 2 reconstructs w on the dual volume centroids in
x and periodic z, using MC slopes and the existing relative positive-trace
limiter. This modifies face reconstructions only; it never floors U or n.
The mesh is not uniform in x and the axis/wall dual supports are half cells.
Impermeable radial flux is required explicitly. Face-positive reconstruction
is not a theorem of positive completed updates; the caller checks both RK
stages and fails atomically if the species partition is invalid.

This is a diffusive transport scheme, not the spatial logarithmic mean of the
fixed-density electron-current split. It conserves integral U, not exact
physical logarithmic entropy. Numerical entropy/energy partition errors must
be measured. `AdiabaticSpeciesMHD` owns the shared two-stage integration.
