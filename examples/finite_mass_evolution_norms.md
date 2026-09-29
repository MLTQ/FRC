# Manufactured trajectory errors and changing signals

Independent physical cell averages of both parent species provide conservative
initial/final references. Electron pressure uses its nodal representation; ion
pressure is the conservative cell partition. Reports absolute RMS/max/axis
errors, final relative error and initial representation error separately.

Error in the actual numerical change is also normalized to the true continuum
change, avoiding a large background total energy or pressure hiding evolution
error. Subtracting the initial representation error is explicit and does not
replace final-state error; both are retained. Quadrature-order replays measure
reference integration uncertainty. This module never changes a numerical state.
