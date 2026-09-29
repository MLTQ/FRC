# Fixed-grid time comparisons

Reads completed archived trajectories at identical mesh and endpoint, compares
physical B/J and Q1 electron pressure by independent quadrature, and cell ion
pressure by its defined representation. Field differences are evaluated directly
to avoid subtracting a large stationary magnetic bias.

Returns dimensional RMS differences and successive halving ratios. Tiny ion
signals may be limited by subtractive thermal recovery; limiter switching and
the fractional last wave step may affect ratios. No fourth-order rate is forced
or inferred merely from classical RK4. This helper never reruns or edits a state.
