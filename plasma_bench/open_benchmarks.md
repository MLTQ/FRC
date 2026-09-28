# Independent asymmetric open-boundary checks

`CurrentPacket` is a C3 compact toroidal current displaced from the midplane;
it excites both odd and even magnetic multipoles. `reference` integrates exact
circular-loop fields over the smooth cross-section using tensor Gauss quadrature.
It includes an exterior axis probe. Chosen probes lie outside current support,
avoiding singular quadrature. This is a known source, not a dynamical equilibrium.

`static_packet` solves Kq=mu0 J dA with the open exterior and returns unknowns
and source vector for independent residual/field checks. `nodal_probes` requires
exact node alignment, avoiding interpolation error in spatial comparisons.
The outside probes use analytic multipole derivatives, not field-grid gradients.
Mesh, modal truncation, quadrature and matching-surface errors are separate.
