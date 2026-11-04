# iec_diagnostics.py

Full live-ion statistics for HOPS and the positive-cage model. Equal macro
weights are inherited from `Particles`. No species, velocity or position is
changed; the controller only receives an observable.

R²=mean(|x-xbar|²), Rdot=mean((x-xbar)·(v-vbar))/R. Principal covariance
lengths, centroid, origin radius and 10/50/90% radii expose departures from a
sphere. The RMS scalar does not impose spherical or axisymmetric evolution.

For a=mean(dx·dv)/mean(dx²), decompose ion kinetic energy exactly into centre
translation, homologous breathing (a dx), and residual motion. The last term
includes shear/rotation and is **not** a thermodynamic ion temperature. Report
the independent sum-minus-total identity without repairing it. Absorption
changes the observed population: radius changes alone do not prove compression.
`ion_sensor` omits eigensolves/quantiles for per-step feedback but uses the same
radius and radial-speed definitions. No empty-population radius is fabricated.
