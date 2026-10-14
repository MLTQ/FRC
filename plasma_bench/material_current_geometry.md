# Moving current geometry with general species temperatures

`meridional_geometry` constructs the same physical Piola maps, particle-weighted
carrier mass, magnetic energy, Faraday coupling and actual metric derivatives
as MaterialCarrierCurrent. It accepts general entropy/temperatures because the
thermal subclass must evolve their pressure work and transport explicitly.
The isentropic base retains its original restrictions and frozen source while
its higher-order trajectory is live. Consolidation into this shared routine
will follow that trajectory's completion; do not edit its live source closure.

There is no bulk rotation or poloidal magnetic field in this sector. Actual
fluid geometry/thermodynamics are evaluated without a replacement entropy,
density floor, vacuum atmosphere or frozen deformation. The return contract is
the existing current geometry dictionary; no new numerical equations are added.
