# Infinite-exterior report

Runs current-source static, modal, matching-surface, pulse mesh and CFL studies.
Saves complete configurations/histories, dependency versions and a fingerprint of
all package modules plus the three measurement/report dependencies. `plot` shows
independent convergence, modal saturation, fixed/open work sensitivity and the
remaining pulse-flow resolution error. `report` computes sensitivity percentages
from actual results and states their denominator (the 4 m result).

The measured full-space response budget includes positive beyond-grid induced
magnetic energy; it excludes prescribed external self/mutual and circuit energy.
The report does not conflate infinite magnetic space with free fluid interfaces,
transport validity, magnetic-shock accuracy or formation/stability validation.
No partial or failed evolution is reported as a completed benchmark.
