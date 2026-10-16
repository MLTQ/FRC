# Complete joint charge/excitation source evolutions

Combines the implemented components in one evolving finite electron reservoir,
using N=8.961429541213399e19 m^-3 and the earlier 3.5 us argon source protocols:
initially all Ar+, 10% Ar+ seed, and 10% seed with a 20 eV/nucleus heat pulse.
All begin at Te=10 eV; heat is a sin-squared 1 us pulse with an independent
analytic integral. Prescribed heat is a source experiment, not coil coupling.

The joint network has 15 neutral levels, 40 Ar II LS terms and ground states for
charges2..18. Electron number, heat and produced photons are evolved alongside
all populations. Neutral and ion ADF04 files plus an Ar II ADF07 file are
explicit personal input paths; no provider tables or reconstructed rate arrays
are archived or downloaded by this example. Input hashes identify the sources.

ADF07 initial IDs1/2/3 map by configuration, multiplicity, L and g to Ar II
2P/4D/4F (local indices0/2/5). Final code0 is unresolved, collapsed to ground
only as the explicit preliminary source approximation. Ground Ar II rates are
replaced, not double counted. Neutral ADF04 S is excluded in favor of the one
existing full-ladder ground DI fit. Other excited ionization and excited-parent
RR are missing and their occupied fractions are saved. RR capture immediately
cascades to ground; it does not populate new resolved levels. No complete
collisional-radiative or optically thick energy closure is claimed.

Joint cases run at700/1400 steps (5/2.5 ns), with matched1400-step controls:
the former ground-only source, and the same joint model excluding the two new
metastable ionization channels. The latter is an explicit physical-channel
ablation, not an equally credible final model. Ground-only also changes the
Ar II ground rate source; differences must not be attributed solely to excitation.

Six complete existing independent limits evaluate the new coupled integrator:
hydrogen DI with scalar time-versus-ionization quadrature, and RR with exact
finite-reservoir power-law solution. These use the existing mathematical
references, not new tests or a reporting framework. All accepted states and
rejections are retained, with number, charge, energy, external heat and missing
rate population observations. No state correction is performed.

An additional40ns neutral resonance decay starts with exactly zero electrons
and electron energy. It must still emit spontaneously; the independent reference
is exp(-A*t), with photons DeltaE*(1-exp(-A*t)). This exercises the corrected
zero-electron path in the full network without creating artificial carriers.
The earlier outputs/atomic_network source archive is retained; revised results
use a fresh directory, never overwriting the completed previous histories.
The intermediate outputs/atomic_network_zero_carrier attempt rejected its first
neutral midpoint: dense matrix factorization introduced negative1e-18-level
populations in a disconnected, unpopulated charge block. The rejected array and
accepted initial state are retained in the output; the rate-containing generator
is preserved only in the private personal-input directory, with recorded hash.
The next source uses
sparse exponential action, preserving this structural zero without clipping.

The source is homogeneous with fixed nuclei and no fluid/current feedback;
no motion, neutral momentum, absorption or real photon escape is evolved.
Physical formation remains a separate full-goal requirement.

Final run55867 EXIT0 in89.663s, source25eb4512063bdfb10f3395c52d8f85b31f04d6d976854dd341de8b7e26448a33;
all19 trajectories, PNG inspected/PDF saved in outputs/atomic_network_sparse.
Maximum joint time-halving population/temperature differences are1.72801e-5
and0.00104980eV. All-argon number/charge/energy residuals are<=2.72e-13,
1.41e-13 and4.90e-12eV; forcing quadrature error is separately larger.
Two excited channels contribute36.57% of ionization events in the initially
all-Ar+ case, with other populated channels still missing. Full measured scope,
old failure, independent convergence and next requirements: docs/ATOMIC_NETWORK.md.

`examples.bsr_atomic_source` reuses `argon_models` and `evolve` for its matched
full-history comparison, replacing only the ground-neutral ionization provider.
The baseline model and this original example are unchanged.
