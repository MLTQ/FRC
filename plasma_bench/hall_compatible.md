# Curl-compatible nonlinear Hall geometry

Reuses the legacy topology, Er/Ez spaces and work-projection operations, while
replacing Etheta by the psi/r shape and the magnetic energy by the explicitly
declared compatible_magnetic convention. Interpolating psi linearly in r²/2,z
and taking its physical curl now commutes with the coefficient Faraday map.
At the regular axis Etheta and Br are proportional to r. No fitted coefficient
or special Hall-rate correction is applied.

The same new Etheta interpolation builds its current and electric mass. The
magnetic gradient uses the fixed compatible energy, independent of electric
quadrature order. `energy_density` and `energy_rate` expose that convention for
pressure recovery and local Poynting verification. Toroidal <r²> g² is unchanged.
Default common Hall quadrature is order 5; controls can vary it independently.
Construction replaces the inherited theta factor once; all application-time
operations use only the compatible maps/factor. The legacy component remains
available to reproduce its diagnosed 21.46% axis defect.

This class is not an implicit change to ToroidalState. Consumers must explicitly
select the compatible energy; energy projection is not a thermal repair.
