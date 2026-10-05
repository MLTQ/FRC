# Wall-conforming continuum fixture initialization

Wraps the frozen analytic contact/wave initializer. Dimensions, physical masses,
density, temperatures, point sampling and poloidal flux are unchanged. The
outer normal-fluid slot is an inactive exact zero. Analytic J1(qR) floating
residue is omitted at initialization. The cell-integrated toroidal approximation
is projected once with the explicitly declared magnetic-mass boundary projector.
Both initial and evolved physical field errors must be measured independently.

Returns the original five-item fixture tuple with StrongWallFluid, geometry and
plasma. `system.original_system` is retained solely to construct an approximate
preconditioner and diagnostic controls; it is never used as the evolved RHS.
This is a bounded, periodic verification fixture, not a lab operating point.
