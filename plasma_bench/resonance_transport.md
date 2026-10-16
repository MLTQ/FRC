# Finite-cylinder resonance photon transport

`transport` follows a weak pulse of resonance quanta, initially excited atoms
at the cylinder centre, until every packet crosses the side or an end. The
uniform stationary ground-neutral density is prescribed only inside the finite
cylinder; exterior opacity is exactly zero. Surface encounters remove photons
from the gas. This is not a quartz transmission/reflection model or a neutral
atmosphere added to a formation run. No weighted killing, escape factor,
frequency cutoff, spatial mesh or transport timestep is used.

For a single ground-resonance branch with rate A, each excited residence time
is Exp(A), followed by an isotropic photon. Flight distance is Exp(kappa_nu)
or the exact ray-cylinder boundary distance, whichever is shorter. Flight time
is distance/c. Absorption returns the quantum to an excited atom at that point.
The low-excitation limit keeps ground density fixed and neglects stimulated
emission, quenching, branching, continuum absorption and ionization. Counts of
excited atoms, in-flight photons and escaped photons are accumulated at supplied
output times by exact interval endpoints; output spacing is not a solver step.
Each packet carries the same h nu0 energy (narrow-line approximation). Frequency
shifts do not exchange kinetic energy here; largest relative detuning is saved.

`redistribute` retains coherence in the atom's rest frame. Set x=(nu-nu0)/sigma,
a=gamma/sigma. The absorbing atom velocity parallel to the incoming ray has
conditional density proportional to exp(-u²/2)/[(x-u)²+a²]. Transverse velocity
components are independent unit normals. The new lab detuning is
`x_new = x_old - u_parallel + u_vector dot direction_new`.
`absorbing_velocity` samples this distribution with exact rejection envelopes:
a Lorentz proposal in the core; a mixture of truncated normal and truncated
Lorentz proposals in the wings, bounded on either side of a positive split.
The split changes efficiency only. No core/wing redistribution switch is made.

An explicit `complete` comparison resamples the full Voigt emission profile
after every absorption. It is a mathematical closure comparison with the same
line shape, not a prediction of a high-collision gas whose broadening is ignored.
The default atom-frame prescription is a scalar coherent-scattering model, not
a complete quantum redistribution calculation. Isotropy/polarization, stationary
atom positions, fixed natural-lifetime waiting even in far wings, constant bath,
no recoil and no collisions are limitations requiring physical assessment.

The initial spontaneous Voigt frequency is sampled as independent Gaussian plus
Lorentz variates without truncation. Formal Voigt wings extend to infinite
detuning; the maximum relative detuning diagnostic exposes narrow-line limits.
All packets must escape; exceeding the generation guard raises rather than
silently discarding survivors. Per-packet exits, counts, residence/flight times
and complete event lists for the first eight packets are returned.

| Dependent | Contract |
|---|---|
| `examples.resonance_transport` | independent packets, explicit seed, SI geometry/density/time, all terminal exits |
| `RadiativeLevels` | source energy, natural/thermal widths and Einstein line strength |

Primary background: Holstein (1947), DOI10.1103/PhysRev.72.1212;
Rajaraman & Kushner (2004), DOI10.1088/0022-3727/37/13/009;
Alves-Pereira et al., https://arxiv.org/abs/physics/0701209. The latter identifies
limitations of replacing atom-frame redistribution by a core/wings heuristic.
This implementation uses the conditional thermal velocity distribution instead.

The completed 20-case study and limits are in `docs/RESONANCE_TRANSPORT.md`.
Complete randomization produced rare relative detunings up to 0.22337, outside
the narrow-line approximation; those packets are preserved in the evidence.
The independent reference checks free flights, not many-scattering experimental
accuracy. Detector-resolved, distributed-source validation is still required.
