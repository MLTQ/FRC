# Sourced electron–neutral elastic force and heat

The missing ground-neutral/electron channel now has a finite-drift Maxwellian
moment implementation and complete homogeneous relaxation evidence. Read
`plasma_bench/elastic_maxwellian.md` for equations, stationary-target energy and
COM-angle assumptions, and `plasma_bench/elastic_atomic.md` for source coupling.
The binary force uses the momentum-transfer cross section, not a total or
effective-inelastic cross section. Recoil and signed species heat come from
the same elastic events; total pair momentum/energy cancel algebraically.

The physical case uses the already authorized private LXCat BSR-500 Ar ELASTIC
table,0–300.0014eV,462points. Source authors O. Zatsarinny and K. Bartschat;
[provider](https://www.lxcat.net/BSR),
[LXCat description of elastic momentum-transfer data](https://nl.lxcat.net/notes/index.php?download=lxcat_team1).
Input SHA `ed851c2f5f9292ff46ff890069a28b3f4c4edba1068ae130c0a3dd6080758f44`.
Raw personal files remain outside Git. Neither private tables nor reconstructed
rate curves are redistributed. Linear interpolation and finite supplied domain
are explicit approximations; unknown cross sections are never extrapolated.

## Complete independent source evolution

Run7161 EXIT0, `outputs/elastic_pair_relaxation`, source
`012788ed4b91c4e561f1c0c5f3bf63388e8b306ae833b3066cdd0e0b5efe8a6d`.
All six100/200/400-step trajectories complete eight initial small-drift relaxation
times, with every state archived. Separate radial/polar velocity integration
uses no production angular kernel. Full8x32 and16x64 reference quadratures agree
within1.844e-9m/s and1.688e-14eV for electron/argon. No temperature/energy repair.

The equal-mass constant-cross-section verification fixture lasts0.698405ms.
Its zero-drift drag agrees with the independent formula to2.22e-16 relative,
and temperatures approach the independently energy-derived3.75eV equilibrium.
Finest drift error2.145e-9 of initial, temperature error5.416e-8eV and total
energy drift1.436e-8; timestep-halving decreases errors about16fold.

The sourced electron/argon case uses physical electron and39.948u masses,
n=1e18/1e20m^-3, initial T=2/0.02585eV and relative drift three combined thermal
speeds. Its complete duration is0.707794548us. There is no field sustaining the
drift, no reactions and no spatial/chamber dynamics.

| RK4 steps | Drift error / initial drift | Maximum Te error, eV | Maximum relative total-energy drift |
|---:|---:|---:|---:|
|100|2.59186e-5|4.02608e-4|1.81696e-5|
|200|1.87306e-6|2.73370e-5|1.07745e-6|
|400|1.19766e-7|1.70782e-6|6.42193e-8|

The finest final Te/Tn are7.996615603/0.02588385442eV; drift is0.340392m/s.
Frictional drift energy mostly heats electrons over this short interval; the
pair has not reached the eventual0.104802eV common equilibrium. Maximum scaled
momentum error is1.05e-12. Maximum probability beyond the supplied table is
3.61e-16; this is not a bound on missing force/heat because the unknown cross
section and extra velocity weights matter. Scientific PNG inspected; PDF saved.

## Moving coupling and remaining work

The explicit `--neutral-elastic-file` option now adds this single ground-neutral
pair to the moving source, alongside every occupied charged pair and the existing
reactions. Both force and heat enter the same conservative pair matrix and
collision-half/reaction-full/collision-half integration. Actual populations and
physical masses remain; neutral momentum is evolved rather than held in a bath.
Full20ns evaluation79419 is terminal EXIT0 in
`outputs/reactive_material_argon_neutral_elastic_10ns`:971steps (527/444),
21349density cells,103.579s,Erel2.84185e-11. Source SHA
`6539c5fcd1de3d50fc76e1ad923341d0a67808b1672f2858c5b2845d4628b8e1`.
The supplied table is used only for occupied electron/ground-neutral pairs.
Maximum evaluated relative-velocity probability outside its domain is9.304e-64.

Existing full-endpoint comparison99437 EXIT0 gives E/J relative differences
1.02593e-8/6.76392e-5 from the matched charged-collision case. Electron front
shift4.24544e-13m, common-volume velocity/T RMS0.00562399m/s/3.01909e-5eV;
neutral velocity/T RMS1.99583e-7m/s/2.11014e-8eV. These norms include displaced
sharp features and do not represent pure integrated collision heating or current
attenuation. The n1e14m^-3 slab remains a numerical component; current/thermal
and source-time accuracy are unresolved. Bead frc-transport.7 remains open.

Missing ion-neutral/neutral-neutral channels, charge exchange, excited-target
scattering, molecular channels, distribution evolution and physical screening
remain. The local Maxwellian closure is a modeling assumption, especially at
finite drift. The sampled phase-space moments do not prove Maxwellian evolution
or validate a laboratory discharge. No unit suite or publication occurred.
