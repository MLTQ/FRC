# Independent Hall integration controls

Adds two lower-amplitude slow Hall waves on the same nr16 mesh and same phase,
separating linear-reference accuracy from finite-amplitude residuals. Uses the
unchanged frozen engine/measurement code and records this additional protocol's
own source hash in every configuration and collected record.

The resistive reference is Btheta=C J1(qr) exp[-(eta/mu0)q²t], J0(qR)=0,
uniform positive argon carrier density and an axial bias. This is a truly
conducting resistive trace: Jz(R)=0. z is uniform, so Hall induction vanishes
in the continuum. The finite toroidal magnetic force creates O(C²) fluid and
poloidal changes; C=1e-7 T makes its feedback on Btheta/heating negligible at
the measured spatial scales, not mathematically absent.

Independent exact annular g averages and continuum magnetic-energy loss give
the expected field and total eta J² heat. The pressure change is not compared
to heat alone: compressional pressure work is also O(C²). All actual coupled
steps, positive pressure, heat and nonlinear residual histories are retained.
Four mesh cases and three extra time caps measure spatial/current/heat and
temporal behavior separately. No electron/ion heat partition is inferred.
