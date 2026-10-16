# Resonance-photon trajectories through the provisional tube

`python -m examples.resonance_transport NEW_OUTPUT --packets 4096` follows two
NIST Ar I resonance lines through a 1 m long, 15 cm diameter neutral cylinder
at 300 K. All initial excitations are at its centre; this weak-pulse Green
function calculation is not a self-consistent discharge, formation experiment,
uniform volume source or specified lab preparation. Exterior opacity is zero.

Two vacuum runs and sixteen absorbing runs span 0.1/1/10/100 mTorr equivalent
ground-neutral density. Each packet is tracked until an actual side/end crossing.
The default uses conditional atomic velocities and coherent atom-frame frequency;
a complete-randomization comparison uses the same Voigt absorption profile.
Two additional 106.7 nm, 10 mTorr runs use one-quarter/four-times as many packets
with separate seeds to measure statistical sensitivity. Means include sample
standard errors; stochastic runs are not claimed to have deterministic dt order.

Independent references integrate Beer-Lambert transmission over emission
frequency and direction, without ray stepping or sampled paths. The vacuum
reference also integrates the exponential emission delay plus exact flight time.
For side radius R and half-length Z, mu*=Z/sqrt(R²+Z²), and the exact mean
boundary path is R asin(mu*)-Z log(mu*). Full escape CDFs and first-flight escape
probabilities are compared, including statistical errors. These checks do not
establish the validity of the redistribution prescription in the laboratory.

Output times only collect exact event intervals; they are not integration steps.
All per-packet escape times, exit positions/frequencies, absorption counts and
total photon-flight times are saved, together with interior excitation/photon
counts and full event lists for the first eight packets. Energy counting uses
constant h nu0 per packet; velocity recoil and frequency-shift energy exchange
are excluded, and largest detuning/nu0 is reported. Conservation of this packet
inventory is narrower than conservation of a coupled plasma/atom/photon system.

The source/input archive and PNG/PDF use the existing project record tools.
No unit tests, escape-factor fitting, tail truncation, artificial gas outside the
cylinder, stochastic killing or photon-budget normalization are introduced.
Pending physical work includes quenching, metastable/ionization channels,
moving nonuniform species, atom transport, polarization and wall response.

Completed run 92249: 20 cases, 524.793 s, all packets escaped. Measured values,
statistical sensitivity and formal-wing limitations are recorded in
`docs/RESONANCE_TRANSPORT.md`. The archived output source remains the exact code
used for this run. Published detector decay rates need matched source/detector
geometry; they must not be compared to inverse mean global escape time.
