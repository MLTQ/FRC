# Continuous full-trajectory energy audit

Independently reconstructs initial/final n,v, both C/n, J and magnetic increments
at physical quadrature orders 16 and 20. Integrates bulk kinetic, both physical
internal, current and magnetic energy changes separately. The full local
constant-bias cross term is retained and measured rather than assumed zero in
the audit. Stable expm1/log1p thermal differences preserve small perturbations
against the large background without any correction to states or energy.
Linear scalar and bias-cross reconstructions act on conservative changes before
background subtraction; ds=(S delta C-s0 S delta n)/n. This exact algebraic
evaluation avoids avoidable interpolation/subtraction roundoff. The initial
direct-background audit is retained in regular_plasma_energy_before_difference.

Reports wave-signal-scaled energy change, particle drift, species entropy changes
and temperature extrema separately from the original grid energy change.
Hashes original raw arrays, freezes the complete coupled closure plus this
independent audit module, and refuses overwrite. This integrated audit does not
prove local geometric enthalpy/Poynting interface conservation.
