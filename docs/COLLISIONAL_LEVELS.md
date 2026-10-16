# Resolved excitation with a finite electron-energy reservoir

Run 61889 completed all 16 cases in 13.710 s, exit 0. Full accepted population,
electron-energy and produced-photon histories, source archive, input hashes and
both inspected PNGs are in `outputs/collisional_levels`. Source:
`4096b6aae2627eef2c7ee1f8ba3d4c99310305bd319f01c64ce1d0b5d22c01db`.
This adds a missing atomic source component; moving-plasma coupling remains open.

## Physical inputs and scope

The [OPEN-ADAS neutral argon input](https://open.adas.ac.uk/element/detail/adf04/arlike/arlike_mom08%5D%5Bar0.dat)
supplies 15 levels and 91 Maxwellian excitation pairs over 500–500000 K. Its
comments attribute digitized, Maxwellian-convolved transitions to Dasgupta et al.,
[PRA 61, 012703](https://doi.org/10.1103/PhysRevA.61.012703) and the
[erratum](https://doi.org/10.1103/PhysRevA.65.039905). The component explicitly
excludes the supplied ground ionization S row. It uses the repository's current
NIST 15-level energies and radiative A values after matching ordered weights
and energies within 1 cm^-1. The older supplied resonance A values differ
substantially; this replacement is explicit, not a claim they are identical.

The [OPEN-ADAS singly ionized argon input](https://open.adas.ac.uk/element/detail/adf04/cllike/cllike_dcg07%5D%5Bar1.dat)
provides 40 LS terms, 780 pairs and 4000–800000 K. It attributes the data to the
452-state RMPS calculation of Griffin and Ballance,
[DOI 10.1088/0953-4075/40/23/013](https://doi.org/10.1088/0953-4075/40/23/013).
Its source xj encodes (g-1)/2, not the physical J of a fine-structure level.
Source term energies and term-averaged radiative rates are retained. No resolved
Ar II line-profile or opacity accuracy follows from this term model.

Both files are explicit personal local inputs in
`/private/tmp/frc-openadas-personal/`; neither raw tables nor reconstructed rate
arrays are bundled in the repository or output archives. See the
[provider terms](https://open.adas.ac.uk/terms-and-conditions) before any future
redistribution or incorporation into a published system. The parser does not
download data. Exact identities and paths are in the saved summary:

- Neutral `arlike_mom08_ar0.dat` SHA256:
  `20f79ac06aa08e191f3a647a42734abca02651e5751248107c6928a52a6da2b4`.
- Ion `cllike_dcg07_ar1.dat` SHA256:
  `c17c12fafd2679c1e229c79a2d0919cdaa1096a0d3491de26ff1d64df04a1af9`.
- NIST radiation SHA256:
  `bf04549cd6dae143ff65ee92700902f1eea0567003b005fea5310dd598d5f2db`.

The homogeneous component has fixed ne=1e20 m^-3 and one electron per modeled
atom. Unevolved background ions balance the neutral component's electrons.
There is no charge change, density evolution, fluid/current feedback, absorption,
escape, finite-drift correction or non-Maxwellian electron distribution. The
data bounds are enforced without extrapolation. Missing processes are not zero.

## Equations and numerical update

For effective collision strength Upsilon(T), statistical weights g and level
spacing DeltaE, with T_K in kelvin and T_e in eV:

    q_up   = C Upsilon / (g_lower sqrt(T_K)) exp(-DeltaE / T_e)
    q_down = C Upsilon / (g_upper sqrt(T_K))
    C = sqrt(8/(pi*m_e))*pi*a0^2*Rydberg/sqrt(kB)
      = 8.629132163236302e-12 in SI.

Upsilon uses linear interpolation in log temperature with nonnegative endpoint
weights. Actual excitation and de-excitation transfers evolve every supplied
pair. Their electron-energy source is -DeltaE*(R_up-R_down). Spontaneous emission
adds A*DeltaE*f_upper to a separate photon-production reservoir. Thus population
sum and electron + atomic + produced-photon energy are generator invariants.
Electron energy is independently advanced, never reconstructed or corrected.

An exponential half-step predicts the midpoint electron temperature, followed
by a full exponential at that temperature. This is second-order exponential
midpoint for nonlinear feedback, not an exact nonlinear evolution. Negative or
nonfinite state, nonpositive electron energy and out-of-range temperature reject;
no clipping, population normalization or energy repair is used. See the module
companions for parser and constructor contracts.

## Completed 3.5 us histories

Each initial temperature ran at 350 and 700 steps (10 and 5 ns). Fine endpoints:

| Component | Initial Te, eV | Final Te, eV | Produced photons, eV per atom |
|---|---:|---:|---:|
| Ar I |5|2.908135|3.119037|
| Ar I |10|3.274114|10.058550|
| Ar I |20|3.775686|24.285135|
| Ar II terms |5|2.206923|4.174986|
| Ar II terms |10|2.249370|11.609106|
| Ar II terms |20|2.288991|26.547388|

The maximum time-halving population L1 difference across these six histories
is 4.28399e-5; maximum temperature difference is 4.31568e-4 eV. Maximum relative
total-energy change is 5.89306e-13 and population-sum error 1.15019e-12 across
both step sizes. Two timesteps measure sensitivity, not an observed order.
Matter entropy decreases while energy transfers into produced photons; photon
entropy is not modeled, so this is not a closed-system entropy measurement.

Produced radiation is not established escape or net laboratory cooling.
Previously measured neutral resonance trapping makes that distinction material.
These fixed-charge curves cannot be substituted for the full argon plasma:
ionization, higher levels, reabsorption and the finite-frequency transfer model
can change the result. The old radiation far-wing limitation remains unresolved.

## Independent closed equilibrium

Spontaneous radiation is explicitly disabled for a separate reference. Each
network evolves 1 ms from ground population and Te=10 eV, using 10000/20000
steps (100/50 ns). The independent thermodynamic solution solves

    15 eV = 1.5*T + sum(E_i*g_i*exp(-E_i/T))/sum(g_i*exp(-E_i/T)).

It uses no collision-matrix eigenvector or fitted rate. Fine Ar I reaches
4.15648524605 eV, within 1.38e-13 eV of the reference; population L1 error is
4.05673e-11. Fine Ar II reaches 4.42304409323 eV, within 1.05e-10 eV;
population L1 error is 6.82963e-11. Across the four runs, maximum energy drift
is 2.37169e-11, population-sum error 6.59246e-11 and smallest sampled matter
entropy increment -8.35e-14 kB, retaining accumulated floating-point effects.

Transient time-halving differences are larger than those endpoint errors:
population L1 maxima 4.56752e-4 / 0.00391361 and temperature maxima
0.00161121 / 0.0242344 eV for Ar I / Ar II. An accurate equilibrium endpoint
does not prove equally accurate transient evolution.

Next couple excitation to a consistent charge/electron-number/energy source,
with consistent ground thresholds and explicit missing excited-state processes.
Then connect photon absorption/transfer and the moving model. Do not add this
as an unqualified optically thin energy sink. Bead frc-transport.5.3, parent .5,
the numerical/physical acceptance gates and the entire FRC goal remain open.
