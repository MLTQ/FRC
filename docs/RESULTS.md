# Initial experimental findings

These are results of a prescribed, stationary conductor model, not predictions of stable plasma formation.

## Reproducible comparison

Chamber radius 5 cm, length 30 cm; plasma radius 3.5 cm, length 20 cm.
Two one-turn filaments: radius 6.1 cm, z=±5.5 cm. Bias +5 mT.
Each waveform has a peak absolute **vacuum center perturbation of 15 mT**.
The unipolar pulse has a 5 µs cosine rise, 10 µs hold and 15 µs cosine fall.
The damped sinusoid is 50 kHz with 20 µs amplitude-decay time; the burst is
three 50 kHz cycles with a Hann envelope. These are prescribed I(t), not circuit models.
Conductivity is fixed at 20 kS/m. Reference material: Ar+, ni=10^20 m^-3,
Te=5 eV, Ti=1 eV, fully pre-ionized. These are illustrative inputs.

Outer domain R=15 cm, |z|≤30 cm, 61×121 nodes; dt=0.25 µs; 100 µs duration.
Topology is sampled every 1 µs. Snapshots do not resolve sub-microsecond lifetimes.

| Waveform | Min Bz/Bbias in plasma | Min on-axis Bz [mT] | Peak abs(Jθ) [MA/m²] | Ohmic heat [J] | Closed-flux lifetime lower bound [µs] |
|---|---:|---:|---:|---:|---:|
| unipolar | 0.955 | 5.000 | 1.453 | 0.1146 | 0.0 |
| rlc | -2.329 | -0.144 | 1.400 | 0.1871 | 1.0 |
| burst | -4.463 | -4.835 | 1.748 | 0.3385 | 10.0 |

- Positive unipolar pulse: no reversal or resolved closed flux for this case.
- Damped pulse and burst: reversal and resolved closed contours occur in the stationary model.
- The oscillatory **vacuum-only fields also reverse**. Reversal alone cannot be attributed to plasma current or called FRC formation.
- The burst exceeds the chosen 5 µs sampled-persistence threshold. This label means only duration of detected closed contours, not stability or confinement.
- Closed contours, O/X points and any axis-connected separatrix estimates are numerical diagnostics at the saved resolution. Ambiguous separatrices and excluded-flux radius remain null.

![Comparison](../outputs/comparison/comparison.png)

![A resolved closed-flux snapshot](../outputs/comparison/burst/closed_flux.png)

## Why these results are not a plasma-formation claim

Reference beta is 9.66, Alfven time 20.21 µs, resistive diffusion time 30.79 µs, and Lundquist number 1.52.
Ion skin depth is 14.34 cm (4.10 plasma radii), and ion gyroradius is 18.20 cm at the bias field.
These ion scales are large; Hall-free magnetized single-fluid dynamics cannot be presumed credible for these parameters.
The run lasts several Alfven times, yet the plasma cannot move. Ohmic heat is comparable to or exceeds the initial thermal-energy estimate, yet temperature and conductivity do not change.
The Drude-inferred electron Hall parameter is only a consistency estimate from arbitrary sigma; it does not supply a validated collision model. Physical collision frequencies and mean free paths are unknown.

## Numerical limits

The independent analytic and manufactured tests pass. The mixed solver's algebraic energy identity closes near floating-point precision, with backward-Euler loss reported separately.
However, temporal and spatial errors are measurable: the weak RLC on-axis reversal changes sign on the coarse mesh. At the nominal mesh, halving dt shifts the minimum axis field significantly relative to that small reversal.
Doubling exterior dimensions at fixed spacing changes the center probe by less than 1% of its peak but affects the near-threshold reversal. The outer boundary is not an exact open boundary.
See [quantified convergence and sensitivity](VERIFICATION.md) before using any amplitude threshold or topology lifetime. The phase diagram is a coarse exploratory map, not a converged boundary.

Joint mesh/time refinement changes the sampled RLC closed-flux lifetime from 1 to 4 µs and the burst from 10 to 5 µs. The existence of short closed-flux episodes survives this check, but their durations are not quantitatively converged.

## Small parameter sweep

27 cases: three waveform families × conductivities [5, 20, 80] kS/m × peak external/bias ratios [1, 3, 6]. Geometry and all other inputs are held fixed.
The four categories separately distinguish no reversal, reversal, closed flux, and closed flux with sampled duration ≥5 µs.

![Phase map](../outputs/sweep/phase_diagram.png)

Equal peak external field does **not** mean equal source energy, magnetic work, or heating. No waveform is declared intrinsically better; the normalization and missing physics matter.

## Next milestone

Implement and verify conservative compressible resistive MHD: mass, momentum, total energy, J×B motion, advection and resistive heating. Require force-balance, wave/compression and energy tests before studying dynamic formation.
Develop transport and Hall applicability alongside it; the present example already warns that a scalar-conductivity single-fluid interpretation may fail.

[Physics and sources](PHYSICS.md) · [Architecture](ARCHITECTURE.md) · [Research queue](RESEARCH_QUESTIONS.md)
