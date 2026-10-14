# A configured moving-plasma experiment

`build(config)` connects tube/plasma dimensions, gas inventory, initial species
temperatures, grid and a biased prescribed waveform to the existing complete
moving current/fluid/species/vacuum equations. This is a simulation input path,
not an additional verification framework or a new physical transport closure.

Neutral fill specifies an **assumed preparation inventory**. With chamber
volume V, neutral temperature Tg and pressure p, Nneutral=pV/(kB Tg). The input
must explicitly assume complete single ionization. Ar, He and Ne give one ion
per initial atom; H2 gives two protons per assumed dissociated molecule. Masses
are declared natural-abundance approximations, not isotope-specific chemistry.
No neutral material remains outside plasma in this model.

The compact initial profile is f=(1-r²/R²)(1-z²/L²), n=nc*f^(3/2),
Te/Ti=Tcenter*f. The actual numerical volume integral sets nc so that the full
ion inventory equals the declared preparation inventory. Both initial species
entropies are analytically constant. No density/temperature floor is introduced.
This assumed localization and heating of the gas is not simulated or established
experimentally. Fill pressure is never mislabeled as the local plasma density.

All supplied loops share one current. Their exact summed vacuum-center field
per ampere converts the specified DC bias plus serialized pulse (in tesla) to
that current. Negative pulses can oppose the bias. Existing Unipolar, SineBurst,
DampedSine and sampled waveforms are available; their initial pulse value must
be zero. The initial induced field/current is zero while the exact source bias
penetrates the plasma. A Helmholtz-spaced pair does not imply uniform field over
the entire long tube. No hardware circuit or power supply is designed.

`TubeMaterialCoilCurrent` stops at any RK stage that reaches the tube radius or
either end; it does not pin the edge, bounce plasma, or invent wall/sheath forces.
The tube is nonmagnetic/nonconducting and otherwise does not affect the vacuum
field. Exterior matter and conductivity remain exactly zero. Coil filaments
must lie radially outside the tube; the inherited source/interface check remains.

`examples.material_experiment` and downstream saved-state displays use the same
builder/config. Metadata records the resolved density, inventory, mass, field
normalization and explicit assumptions. Physical transport, preparation,
ionization and formation conclusions require further work.

Independent drives use `source_groups` instead of `source`: each named group has
its own fixed coils, DC centre-field contribution and serialized pulse. Every
group normalizes its own current by its exact vacuum-centre field per ampere;
the physical applied field is the sum over groups at every actual RK stage.
This supports a Helmholtz-like bias pair separately from a longer formation array,
including different timing and polarity. Coils remain abstract prescribed sources.
`source_configs` exposes the declared groups to physical plots without merging
their waveforms. Metadata retains each group and normalization; the legacy scalar
normalization is null when there is no single common current. All initial pulses
must be zero; this does not simulate the preparation of the initial plasma.

Optional top-level integrator="split" selects the existing symmetric implicit
entropy/current diffusion step; omitted or "rk4" retains the original coupled
RK4. Tube contact is checked on all stages for either solver. This choice changes
time integration, not plasma/source inputs or the physical equations. Strong
source behavior remains to be measured; weak-pulse success is not a guarantee.
