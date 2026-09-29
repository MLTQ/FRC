# Independent principal-symbol derivative

Evaluates the physical primitive mass, ion momentum, separate adiabatic
pressures and Faraday curl of the full ideal/Hall/pressure electric field.
It does not import the reduced symbolic matrix. Perturbations are plane sines
whose values vanish at the target point: only their derivative contributions
survive. Complex-step amplitude avoids subtracting large background rates.
Fourth-order spatial differences evaluate gradients and the nested physical
electric curl, retaining actual n/B gradients and cylindrical current curvature.

Compares the five normalized columns with the independently reduced principal
matrix. Wavenumber and spacing are explicit and can be refined independently.
This is a local equation check, not a global evolution or boundary test.
