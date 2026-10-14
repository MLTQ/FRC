# A finite coil pulse with a freely moving heated plasma

Drives the full moving plasma/current/species/vacuum equations from zero induced
field and current. Two prescribed circular loops at r70mm,z=±35mm share a smooth
unipolar waveform: rise1.5us, hold1us, fall2us. Their common amplitude gives10uT
vacuum-center peak by default. These are abstract source inputs, not hardware
specifications or laboratory recommendations. There is no initial DC bias yet.

Plasma starts with the same compact free-cylinder density/temperature fixture as
the material trajectory, R35mm,Lhalf70mm,nc1e18/m³,physical Ar+/electron masses.
All sides/ends move from calculated forces. There is no exterior gas or
conductivity. The exterior follows actual boundary deformation at every stage.

The complete7.121769us run includes the entire pulse and a source-off interval.
Two energy budgets use independently accumulated waveform/linkage work:
plasma+current+induced magnetic energy against -Idot*lambda, and that energy
plus mutual I*lambda against I*lambdadot. Trapezoidal work observations are not
fed back into state integration. All actual fields, heat, current, motion,
source/response quantities and physical residuals are retained, with every
accepted state and frozen source. No unit tests or reporting framework.

Optional `--integrator split` uses SplitMaterialCoilCurrent's symmetric implicit
upstream-diffusion step around the central physical RK4 evolution. It preserves
every accepted state in the same format and records the actual nonlinear solve
and thermal/field exchange observations. The source-response work integral is
null for this method because the pointwise zero-current trace is ambiguous;
direct waveform work -Idot*lambda remains independently integrated. Defaults and
physical parameters are unchanged. This option is evaluated with complete pulses,
not unit tests. See `plasma_bench/material_split_current.md` for equations.

This first coarse trajectory establishes execution/coupling only. Its field,
current and force accuracy need independent physical/time/spatial comparisons.
The inherited upwind method has already shown current time sensitivity. Constant
resistance and conditional exchange are not validated ionization/transport, and
the cold edge is not certified. No reversal, topology, stability, confinement or
formation claim follows from a complete pulse or a small integrated work error.
