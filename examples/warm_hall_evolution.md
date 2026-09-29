# Coupled warm Hall measurements

Advances the full five-stage composition under explicit case fluid CFL and
phase/time caps. Fixed-count controls take exactly that number of steps;
other cases recompute the changing fluid limit. Every accepted time, total and
species energy, entropy, positive minimum, divergence and nonlinear residual
is stored in compact columns. No trial failure is converted to an accepted run.

Wave initial/final errors use independent physical oracles, including warm
species amplitude and phase. Nonlinear trajectories have no claimed exact time
solution and require spatial/time comparisons. Lossless initial/final fluid,
electron and volume arrays allow independent replay. No physical Ohmic heat
or collisional partition is introduced in these eta=0 coupled cases.
