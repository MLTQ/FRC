# Local warm energy and boundary measurements

Pressure cases compare conservative energy divergence, each physical face power
and local ion-energy defect with independent analytic current/thermodynamic
rates. Exact-zero barotropes report null relative-to-zero norms and use the
nonzero physical pressure-work scale explicitly. Axis errors and Gauss12/18
reference uncertainty are retained. Both global budgets and local errors matter.

Zero-trace diffusion uses a true J1-zero Btheta eigenfunction, annular g means
and implicit evolution to exponent .02. Records retain **both** instantaneous
RHS and evolved-state errors, including the nonconverging first-annulus rate.
This is a standalone magnetic boundary test, not a warm resistive/heating
closure. Positive heat, conservative energy flux and a fine two-half-step
control do not erase the derivative defect. Arrays are lossless in SI.
The magnetic solver's tolerance and iteration budget are read from each case.
