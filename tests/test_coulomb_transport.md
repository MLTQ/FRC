# Collision and local-regime checks

Checks SI implementation against independently published rounded NRL cgs
collision/resistivity and electron inertial-length coefficients. Explicit
Coulomb logarithm numerical cases cover both temperature branches. Known
carrier drift and mass-ratio skin lengths verify local unit conventions.

Unavailable classical logarithms remain null with a reason, not floored;
invalid vacuum/thermal/charge inputs fail. Stage 0/1 reference diagnostics
include conditional estimates without modifying prescribed conductivity or
claiming full physical validity. JSON output excludes nonfinite sentinels.
