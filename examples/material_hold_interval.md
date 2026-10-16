# Saved strong-pulse constant-current continuation

Localizes the measured full-system energy loss in an actual source-driven
trajectory. The default source states are steps1450..1500 of the4000-step
energy-diffusion run:1.26875..1.3125us. Both coil currents are exactly constant;
the independent physical limit is constant total fluid+current+induced-field
energy. Source self-energy is excluded under the existing energy convention.

Starts each full continuation from the same saved arrays, without projection,
repair or changing density/temperature/fields. Uses the existing diffusion
halves and full central RK4 functions at50/100/200 steps over43.75ns. Observes
energy before/after each actual operator, retaining all accepted full states and
failures. Only time resolution changes. This local interval cannot replace
whole-pulse time/space convergence or physical transport/edge validity.

The coarse endpoint is compared with its immutable archived endpoint to verify
that this split observation follows the original trajectory. Physical endpoint
B/J/velocity/position comparisons use the existing material_samples evaluator,
10-point Gauss quadrature per reference cell and actual finest-run volume.

At the start, middle and end of the source interval, centered directional
differences of total energy along the actual central derivative are evaluated
at four decreasing time-sized perturbations. These are independent of the
coded instantaneous power-pairing identity. The diffusion's separate energy
increments are observed directly. Neither diagnostic changes a state or supplies
a compensating source. No new unit suite, reporting framework or tolerance gate.

Output freezes this example's imported source and copies exact initial/final
input arrays with hashes. It may read an already-accepted prefix of a still-live
run, but never restarts or modifies that run. Keep any live moving source closure
frozen. Use a fresh destination and preserve earlier failures and full pulses.

Completed run72041 (386.734 s) is in outputs/material_hold_interval. All three
continuations succeed. Energy losses at0.875/0.4375/0.21875ns are
0.116611767/0.006434984/0.000243639mJ; directly observed mixing contributions
are only1.17e-14/2.57e-14/4.00e-14J. Central RK4 produces the local loss.
Compared with the finest continuation, endpoint relative B/J differences fall
from6.0401%/33.1348% to0.190964%/1.016104%. Read MATERIAL_ENERGY_DIFFUSION.md
for assumptions and limits. This diagnosis does not establish full-pulse accuracy.

Optional `--central-integrator dop853_fixed --factors 1 2` compares the same
actual-state interval using the higher-order central advance at the selected
outer timesteps. The original saved source config remains in its source run;
the new output config and summary explicitly record the changed central method.
`central_step` dispatches the selected method. New energy histories use the
method-neutral `central_change_j`; earlier frozen RK4 archives retain their
original `central_rk4_change_j` key. Endpoint differences name the actual finest
factor, never assume it is four. Changing the method is not a restart or repair.

Higher-order run94170 completed both selected factors in296.846s, source
df92f8cb55d196b0e89b3fad787044c4027be7585411ef99ee39c245bdfaab92.
Losses at0.875/0.4375ns are2.79612e-8/4.56020e-10J. Endpoint relative B/J
differences are1.42660e-5/6.30950e-5, and position mismatch7.69328e-9m.
Outputs/material_hold_interval_dop853 in the main repository preserves every
accepted state and inspected figures. This is a local method result; complete
strong/weak pulses remain under study with the same physical inputs.
