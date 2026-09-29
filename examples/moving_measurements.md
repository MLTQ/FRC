# Moving-plasma evolution measurements

Runs a declared 26-case protocol: gas and FRC mesh refinement, a driven uniform
bias, timestep controls and FRC matching/mode/quadrature controls. Every case
starts independently from its stated continuum data. Ten equal diagnostic
intervals are sampled; initial, midpoint and final full fields are archived.

Partial records identify accepted time/steps and remain explicitly incomplete.
Failure records preserve the accepted state rather than forcing success. Reuse
requires a completed record, exact full source/environment fingerprint and
matching configuration. The source is rechecked after each run. Lossless NPZ
archives are hashed. No historical fixed-boundary record is modified.

Time caps are upper bounds: the actual material-wave/motion CFL may be smaller.
Reports must inspect accepted steps before interpreting a nominal time study.
No analytic solution is assumed for driven FRC motion. Equilibrium initialization
does not demonstrate formation, lifetime, stability or confinement.
