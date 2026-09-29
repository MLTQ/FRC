# Moving-material study assembly

Requires all 26 evolution cases and 27 static evaluations complete and matching
their frozen full source/environment/configuration and field hashes. Produces
an ordinary lossless JSON aggregate, a derived scientific Markdown report and
two numerical-field/error figures. Report, analysis, formatter and plot sources
are hashed separately from the frozen solver/measurement records. The renderer's
Matplotlib version is recorded alongside the measured NumPy/SciPy/Python versions.

The report distinguishes continuum errors, discrete work identities, nodal
versus interpolated flux, artificial equilibrium motion, physical source work,
magnetic matching effects and remaining model limitations. Time caps are not
silently treated as actual constant steps. It explicitly quantifies the unresolved
driven-flow mesh difference, and makes no formation, confinement, transport or
long-time stability claim. Visual inspection is a separate required
step after figure generation; incomplete studies cannot be published by this
assembler as completed verification.
