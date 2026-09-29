# Compatible Hall/resistive face energy transport

Computes E cross B/mu0 from common electric coefficients and the declared
curl-compatible magnetic reconstruction. Radial-face products are integrated
exactly in z; axial-face Etheta Br uses the same rational radial shape mass
as the magnetic energy. The regular Er Btheta product uses exact <r²>.
Centered magnetic face traces use physical annular r² centroids.

`energy_rate` returns minus the conservative face-flux divergence. Conducting
radial fluxes vanish, periodic axial fluxes telescope, and no boundary power
is silently added. Local difference from magnetic energy change is a numerical
transfer error to be measured, not repaired or converted to physical Joule heat.
The midpoint integrator uses this flux at midpoint fields/electric coefficients.
