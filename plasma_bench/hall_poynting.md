# Hall Poynting transport audit

Evaluates explicit face-averaged E cross B/mu0 from HallCurl's projected
electric coefficients. Radial faces integrate the product of linear axial
Etheta/Bz profiles exactly, and use Ez Btheta. Axial faces integrate the
linear radial Etheta/Br product with physical r weighting, and <r²> (Er/r) g.
Magnetic fields use centered interpolation to the other face locations;
radial interpolation uses exact annular r² centroids. Conducting radial energy
flux is zero and periodic axial fluxes agree, so its FV divergence conserves
the total energy globally. Axis is regular and has zero radial flux area.

`local_energy_audit` independently computes the existing cell magnetic-energy
rate and minus the Poynting divergence. Their difference is the numerical
local thermal defect a future conservative Hall-only update would introduce.
No difference is corrected, redistributed or installed as physical heating.
Global work cancellation does not make this local defect zero. The component
remains experimental until this defect, induction and boundary errors converge
appropriately. It is not enabled in the compressible solver.
