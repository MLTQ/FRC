# Independent warm field audit

Replays quarter-cell density, Q1 electron energy, ion energy and directly
integrated physical magnetic energy using the earlier independent quarter-cell
loops and physical B integration. Dual-volume comparison allows 8 nr machine
eps relative roundoff from subtracting squared neighboring radii; energy replay
has its own 4e-12 relative criterion. No measured physical error is thresholded.
Checks both endpoint pressure archives and explicit periodic/axis/wall flux
traces. Higher-order analytic conservative averages check wave initialization.

Physical wave norms integrate in x=r²/2,z using independent scalar samples and
Gauss orders 8/12, separately from electric-space order 5. The fluid pressure
reconstruction is the declared production representation; its accuracy is
judged against the independent continuum fields, not assumed exact. All four
species pressure/temperature errors, quadrature sensitivity and reconstructed
ion-pressure minimum are retained. Norms are relative to wave perturbations,
with Pa/eV RMS values. No trajectory or state is changed by this audit.
