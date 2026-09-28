# Infinite magnetic exterior regressions

Closed-form dipole flux and analytic derivatives check normalization, signs and
the axis limit. An independent polar-coordinate volume integral checks exterior
dipole energy all the way to infinity; quadrature refinement and eigenvalues
check the multimode energy matrix. Finite differences verify higher multipoles'
vacuum equation and exact radial decay law.

An asymmetric compact current excites odd/even modes. Independent loop
quadrature checks mesh and mode refinement separately, including fields beyond
the matching rectangle and on the exterior axis. Tests verify exact assembled
matrix symmetry, modal pack/expand and the full vacuum solve residual.
Pack/expand is checked in a norm scaled by the whole field, since tiny high-order
amplitudes can have large relative roundoff while remaining below field precision.

A complete pulse must conserve gas plus full-space induced energy while being
insensitive to moving the matching surface. The energy outside the grid must be
positive and counted. Invalid source-enclosing geometry, trace data and exterior
evaluation positions are rejected. These tests do not establish a free fluid
surface, physical transport or magnetic-shock accuracy.
