# Finite argon radiative cascades and neutral absorption

The spontaneous-radiation source now resolves ground + four 4s + ten 4p
neutral-Ar levels and 31 NIST ASD transitions, instead of sending every excited
atom immediately to ground. The new independent `RadiativeLevels` component
does not yet replace the ground-cascade assumption in the DI/RR source or enter
the moving pulse. Its photon output is production, not a declaration of escape.
Read `plasma_bench/line_radiation.md` for the population/energy equations and
Einstein absorption coefficient, and the atomic data README for provenance.

## Complete trajectories

Run61232 exited0 in2.3004s, `outputs/line_radiation`. Frozen source:
`d8ab123e5bb476a0dded1c9ddd23926ba3f667714af0b265fd45772f9498a2e1`.
Atomic record:
`bf04549cd6dae143ff65ee92700902f1eea0567003b005fea5310dd598d5f2db`.
All accepted populations and individual line photon inventories are archived;
PNG inspected, PDF saved. No unit tests, floors, clipping or energy repair.

The actual 2p5→1s4→ground chain has the independent Bateman solution for
rates4e7 and1.32e8s^-1. Over50ns, maximum population errors at1/.5/.25ns
are4.48882e-7,2.65902e-8,1.61660e-9 (ratios16.88/16.45). This verifies
the finite decay evolution, not collision physics or laboratory behavior.

Three complete3.5us cascades start90%ground and1%in each4p level, an explicitly
prepared component study with no electron collisions. The final metastable
fraction is0.04211054, with0.48806210eV/neutral atom still stored and
0.83261963eV/neutral atom produced as photons, from1.32068174eV initially.
The finest step halving changes populations at common times by<=3.8541e-8 and
emitted energy by<=4.55785e-7eV/atom. Number residual<=6.44e-15 and energy
residual<=1.13e-14eV/atom. The metastables' eventual destruction is absent,
not physically predicted to be zero.

## Neutral columns rule out blanket optically thin line cooling

Natural plus thermal Voigt profiles and the Einstein relations give ground
absorption cross sections at line centre and300K of1.15133e-16m² at104.822nm
and3.04201e-17m² at106.666nm. The corresponding vacuum spontaneous lifetimes
are1.8797ns and7.5758ns. These are derived from NIST probabilities, not separately
measured lifetimes in this study.

For a uniform **ground-neutral** gas at300K and a7.5cm centre-to-side path:

| Neutral pressure equivalent | 104.822nm optical depth | 106.666nm optical depth |
|---|---:|---:|
|0.1mTorr|27.79|7.344|
|100mTorr|27794.6|7343.8|

At the earlier source-study nucleus density8.96143e19m^-3, only1%ground-state
neutrals at300K already give radial depths7.738/2.045 and axial half-tube
depths51.59/13.63. At1eV atom temperature the same1%column gives1.259/.330
radially and8.393/2.198 axially. The state population and atom temperature
therefore matter; fill pressure is not interchangeable with evolved neutral
density or total plasma pressure.

These are line-centre column calculations, not integrated escape fractions.
Photons in the wings can escape, and reabsorption does not by itself determine
the final cooling rate. No escape factor has been inserted into A. Pressure,
Stark/Zeeman broadening, motion, redistribution, wall response and the evolving
spatial populations are absent. These results constrain a proposed optically
thin closure; they do not invalidate numerical conservation of the earlier
explicitly thin RR component, or quantify its continuum absorption.

## Next physical work

Add resolved electron-impact excitation/de-excitation, stepwise ionization,
metastable loss, recombination branching and photon transfer before combining
this with the finite-mass/current solver. Avoid a high-charge substitute for
neutral/Ar II physics. Public CHIANTI's argon list begins at Ar III; an INPTDAT
argon rate set was calculated from steady-field non-Maxwellian distributions
and cannot be substituted directly into the Maxwellian model.

The BSR resolved cross sections were located through LXCat2's official download
workflow. Access stops at explicit terms acceptance (attribution, no third-party
redistribution, no commercial use). The user confirmation request is pending;
the terms were not accepted and the data were not downloaded. This does not
block independent simulation work. Keep the full FRC goal and frc-transport.5
open, and collect the live pulse convergence comparisons when they finish.

Sources: [NIST ASD](https://doi.org/10.18434/T4W30F),
[NIST lines query](https://physics.nist.gov/PhysRefData/ASD/lines_form.html),
[LXCat terms](https://us.lxcat.net/instructions/how_reference.php),
[CHIANTI ion list](https://db.chiantidatabase.org/).
