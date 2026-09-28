# Standing magnetic-shock measurements

Runs both actual 2-D kernels against the independently integrated stationary
cylindrical branches, preserving complete parameters, boundary choices, cell
profiles, diagnostics and wall-clock measurements. Errors use exact annular
overlap weights in r=0.44–0.66 m, z=±0.1 m. Density and B are conserved volume
averages; the pressure reference is recovered from independently integrated gas
energy/momentum so subcell kinetic variance is handled consistently.

The cases refine space, time and declared diffusion coefficient, compare the
retained direct implicit solve, and vary the remote inner/outer/axial extensions.
Shorter wall-control comparisons avoid conflating an unrelated remote vacuum
rarefaction failure with local shock accuracy. Failures are retained with their
accepted state and time, never dropped from the report. Individual case JSONs
are written as each finishes. Explicit reuse requires matching package/measurement
source and environment provenance; otherwise cases rerun. Each case's numerical
heat is kept distinct from zero physical Ohmic heat. Global budgets and local
profile accuracy are separate measurements.
