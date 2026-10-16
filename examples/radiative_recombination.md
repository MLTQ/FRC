# Complete recombining source evolutions

Run `python -m examples.radiative_recombination OUTPUT` with a new destination
and noninteractive matplotlib backend. This is a physical component study,
not a unit suite, moving-plasma simulation or formation/preparation prediction.
It archives source, atomic inputs, every accepted state, PNG/PDF and results.

`evolve` integrates charge fractions, thermal energy and cumulative escaping
photon energy with unmodified RK4. No clipping, normalization or floors. Energy
balance includes atomic binding inventory and independent analytic heat input.

`power_law_reference` is an independent exact solution of RR-only evolution:
sigma(E) ∝ E^(-s-1/2) implies alpha(T) ∝ T^-s and mean captured energy
(3/2-s)T. With x the electron fraction and N fixed, eliminating time gives
T=T0*(x/x0)^(-2s/3). Writing beta=2s²/3 and alpha0=alpha(T0),

    x=x0*[1+(1+beta)*N*alpha0*x0*t]^(-1/(1+beta)).

The escaping energy follows independently from initial minus current thermal
plus binding inventory. Full 30 ms runs at 50/100/200 steps use N=1e20 m^-3,
x0=.8, T0=2 eV, alpha0=1e-18 m³/s, s=.7, chi=13.6 eV. These deliberately
specified power-law coefficients are an analytic fixture, not hydrogen data.

Four argon cases use fixed N=8.961429541213399e19 m^-3. Three retain the prior
3.5 µs, initial 10 eV, 100%/10% seed and prescribed 20 eV/nucleus heat cases.
Their 400/800/1600 complete trajectories isolate the added RR source. A fourth
initially singly ionized 1 eV case evolves for 1 ms at 800/1600/3200 steps to
expose recombination over a longer interval. Ionization-only controls are
computed with identical inputs. Cold-source evolution still omits excitation,
DR, three-body processes, molecular channels, transport and opacity; it is not
a realistic afterglow prediction.

The downloaded author's original total RR rates supply an external comparison
for 228 rounded temperature/rate pairs across 12 stages. Four other modern
stages have no original-rate rows in the downloaded master and are explicitly
listed as missing. Rate-fit derivative accuracy is not inferred from this
rate-value comparison. Unknown legacy fit ranges and known extrapolation
content are reported separately; neither is hidden by excluding populations.
